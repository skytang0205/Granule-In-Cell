#include "Pressure.h"

#include "BiLerp.h"

#include <amgcl/amg.hpp>
#include <amgcl/backend/eigen.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>

namespace Pivot {
    Pressure::Pressure(StaggeredGrid const & sgrid):
        m_Grid2Mat(sgrid.GetCellGrid()), m_Pressure(sgrid.GetCellGrid()), m_GradPressure1(sgrid.GetFaceGrids()),
        m_GradPressure2(sgrid.GetFaceGrids()), m_DeltaPos(sgrid.GetFaceGrids()) {}

    void Pressure::Project(
        SGridData<double> &      velocity,
        GridData<double> const & levelSet,
        Collider const &         collider,
        GridData<double> const & fraction) {
        m_GradPressure1.SetZero();
        SetUnKnowns(levelSet);
        if (m_Mat2Grid.empty()) return;
        BuildProjectionMatrix(velocity, levelSet, collider, fraction);
        SolveLinearSystem();
        ApplyProjection(velocity, levelSet, collider);
    }

    void Pressure::Correct(
        std::vector<Particle> &  particles,
        GridData<double> const & levelSet,
        Collider const &         collider,
        GridData<double> const & fraction,
        GridData<double> const & density) {
        m_GradPressure2.SetZero();
        SetUnKnowns(levelSet);
        if (m_Mat2Grid.empty()) return;
        BuildCorrectionMatrix(density, collider, fraction);
        SolveLinearSystem();
        ApplyCorrection(particles, collider);
    }

    void Pressure::InitRestDensity(int numPartPerCell, std::vector<Particle> const & particles, GridData<double> & density) {
        m_NumPartPerCell = numPartPerCell;
        for (auto const & particle : particles) {
            for (auto const [cell, weight] : BiLerp::GetWtPoints(density.GetGrid(), particle.Position)) {
                density[density.GetGrid().Clamp(cell)] += weight / m_NumPartPerCell;
            }
        }
    }

    void Pressure::BuildProjectionMatrix(
        SGridData<double> const & velocity,
        GridData<double> const &  levelSet,
        Collider const &          collider,
        GridData<double> const &  fraction) {
        std::vector<Triplet<double>> elements;

        for (int r = 0; r < m_RdP.size(); r++) {
            Vector2i const cell      = levelSet.GetGrid().CoordOf(m_Mat2Grid[r]);
            double         diagCoeff = 0;
            double         div       = 0;
            for (int i = 0; i < Grid::GetNumNeighbors(); i++) {
                Vector2i const nbCell   = Grid::NeighborOf(cell, i);
                auto const [axis, face] = StaggeredGrid::FaceOfCell(cell, i);
                int const    side       = StaggeredGrid::FaceSideOfCell(i);
                double const weight     = 1 - collider.GetFraction()[axis][face];
                if (weight > 0) {
                    int const c = m_Grid2Mat[nbCell];
                    if (c >= 0) {
                        diagCoeff += weight * (1 + fraction[nbCell] / fraction[cell]) * 0.5;
                        elements.push_back(Triplet<double>(r, c, -weight * (1 + fraction[nbCell] / fraction[cell]) * 0.5));
                    } else {
                        double const theta    = levelSet[cell] / (levelSet[cell] - levelSet[nbCell]);
                        double const intfCoef = 1. / std::max(theta, .001);
                        diagCoeff += weight * intfCoef * (1 + fraction[nbCell] / fraction[cell]) * 0.5;
                    }
                    div -= side * weight * velocity[axis][face] * (1 + fraction[nbCell] / fraction[cell]) * 0.5;
                }
                if (weight < 1) { div -= side * (1 - weight) * collider.Velocity[axis][face]; }
            }
            m_Rhs[r] = div;
            elements.push_back(Triplet<double>(r, r, diagCoeff ? diagCoeff : 1.));
        }

        m_MatL.setFromTriplets(elements.begin(), elements.end());
    }

    void Pressure::BuildCorrectionMatrix(
        GridData<double> const & density, Collider const & collider, GridData<double> const & fraction) {
        std::vector<Triplet<double>> elements;

        for (int r = 0; r < m_RdP.size(); r++) {
            Vector2i const cell       = density.GetGrid().CoordOf(m_Mat2Grid[r]);
            double         diagCoeff  = 0;
            bool           isBoundary = false;
            for (int i = 0; i < Grid::GetNumNeighbors(); i++) {
                Vector2i const nbCell   = Grid::NeighborOf(cell, i);
                auto const [axis, face] = StaggeredGrid::FaceOfCell(cell, i);
                int const    side       = StaggeredGrid::FaceSideOfCell(i);
                double const weight     = 1 - collider.GetFraction()[axis][face];
                if (weight > 0) {
                    int const c = m_Grid2Mat[nbCell];
                    diagCoeff += weight * (1 + fraction[nbCell] / fraction[cell]) * 0.5;
                    if (c >= 0) {
                        elements.push_back(Triplet<double>(r, c, -weight * (1 + fraction[nbCell] / fraction[cell]) * 0.5));
                    } else {
                        isBoundary = true;
                    }
                }
            }
            m_Rhs[r] = 1. - (isBoundary ? 1. : std::clamp(density[cell] / fraction[cell], .5, 1.5));
            elements.push_back(Triplet<double>(r, r, diagCoeff ? diagCoeff : 1.));
        }

        m_MatL.setFromTriplets(elements.begin(), elements.end());
    }

    void Pressure::SetUnKnowns(GridData<double> const & levelSet) {
        m_Grid2Mat.SetConstant(-1);
        m_Mat2Grid.clear();

        int n = 0;

        ForEach(levelSet.GetGrid(), [&](Vector2i const & cell) {
            if (levelSet[cell] <= 0) {
                m_Grid2Mat[cell] = n++;
                m_Mat2Grid.push_back(levelSet.GetGrid().IndexOf(cell));
            }
        });

        m_MatL.resize(n, n);
        m_RdP.resize(n);
        m_Rhs.resize(n);

        m_RdP.setZero();
    }

    void Pressure::SolveLinearSystem() {
        using Solver = amgcl::make_solver<
            amgcl::amg<amgcl::backend::eigen<double>, amgcl::coarsening::smoothed_aggregation, amgcl::relaxation::spai0>,
            amgcl::solver::bicgstab<amgcl::backend::eigen<double>>>;
        Solver solve(m_MatL);
        auto const [iters, error] = solve(m_Rhs, m_RdP);
        std::cout << fmt::format("{:>6} iters", iters);
    }

    void Pressure::ApplyProjection(SGridData<double> & velocity, GridData<double> const & levelSet, Collider const & collider) {
        ParallelForEach(velocity.GetGrids(), [&](int axis, Vector2i const & face) {
            Vector2i const cell0  = StaggeredGrid::AdjCellOfFace(axis, face, 0);
            Vector2i const cell1  = StaggeredGrid::AdjCellOfFace(axis, face, 1);
            double const   weight = 1. - collider.GetFraction()[axis][face];
            if (weight == 0.) return;

            int const id0 = m_Grid2Mat[cell0];
            int const id1 = m_Grid2Mat[cell1];

            // double const p0 = id0 >= 0 ? m_RdP[id0] : 0;
            // double const p1 = id1 >= 0 ? m_RdP[id1] : 0;

            if (id0 < 0 && id1 < 0) return;

            if (id0 >= 0 && id1 >= 0) {
                velocity[axis][face] -= (m_RdP[id1] - m_RdP[id0]) * weight;
                m_GradPressure1[axis][face] = (m_RdP[id1] - m_RdP[id0]) * weight;
            } else {
                double const phi0  = levelSet[cell0];
                double const phi1  = levelSet[cell1];
                double const theta = phi0 / (phi0 - phi1);

                double const pInner = id0 >= 0 ? m_RdP[id0] : m_RdP[id1];
                double const pIntf  = 0.;

                double const intfCoef = 1. / std::max((id0 >= 0 ? theta : 1 - theta), .001);

                velocity[axis][face] -= (id0 >= 0 ? +1 : -1) * (pIntf - pInner) * intfCoef * weight;
            }
        });
    }

    void Pressure::ApplyCorrection(std::vector<Particle> & particles, Collider const & collider) {
        ParallelForEach(m_DeltaPos.GetGrids(), [&](int axis, Vector2i const & face) {
            Vector2i const cell0  = StaggeredGrid::AdjCellOfFace(axis, face, 0);
            Vector2i const cell1  = StaggeredGrid::AdjCellOfFace(axis, face, 1);
            double const   weight = 1. - collider.GetFraction()[axis][face];
            if (weight == 0.) {
                m_DeltaPos[axis][face] = 0;
                return;
            }

            int const id0 = m_Grid2Mat[cell0];
            int const id1 = m_Grid2Mat[cell1];
            if (id0 < 0 && id1 < 0) {
                m_DeltaPos[axis][face] = 0;
                return;
            }

            double const p0 = id0 >= 0 ? m_RdP[id0] : 0;
            double const p1 = id1 >= 0 ? m_RdP[id1] : 0;

            m_DeltaPos[axis][face]      = (p1 - p0) * weight * m_Pressure.GetGrid().GetSpacing();
            m_GradPressure2[axis][face] = (p1 - p0) * weight * m_Pressure.GetGrid().GetSpacing();
        });

        static std::random_device         rd;
        static std::mt19937               gen(rd());
        static std::normal_distribution<> dis(0.0, .01);

        tbb::parallel_for_each(particles.begin(), particles.end(), [&](Particle & particle) {
            particle.Position += BiLerp::Interpolate(m_DeltaPos, particle.Position)
                + Vector2d(dis(gen), dis(gen)) * m_DeltaPos.GetGrids()[0].GetSpacing() * 0.2;
        });
    }
} // namespace Pivot
