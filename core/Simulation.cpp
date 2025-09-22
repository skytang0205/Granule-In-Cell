#include "Simulation.h"

#include <algorithm>
#include <mutex>
#include <random>

#include "Advection.h"
#include "BiBSpline.h"
#include "BiLerp.h"
#include "CSG.h"
#include "Extrapolation.h"
#include "Redistancing.h"

#include <chrono>

namespace Pivot
{
    Simulation::Simulation(StaggeredGrid const &sgrid, double particleRadius) : m_SGrid{sgrid}, m_Collider(m_SGrid), m_Pressure(m_SGrid), m_Velocity(m_SGrid.GetFaceGrids()),
                                                                                m_ConvectVelocity(m_SGrid.GetFaceGrids()), m_LevelSet(m_SGrid.GetCellGrid(), std::numeric_limits<double>::infinity()),
                                                                                m_VelDiff(m_SGrid.GetFaceGrids()), m_Density(m_SGrid.GetCellGrid()), m_RestDensity(m_SGrid.GetCellGrid()),
                                                                                m_DEMGrid(m_SGrid.GetCellGrid()), m_CouplingForce(m_SGrid.GetFaceGrids()), m_TargetFraction(m_SGrid.GetCellGrid()),
                                                                                m_DEMForce(particleRadius), m_ParticleRadius(particleRadius), m_SupportRadius(2 * particleRadius),
                                                                                m_ParticleVolume(std::numbers::pi * particleRadius * particleRadius) {}

    void Simulation::Describe(YAML::Node &root) const
    {
        root["Dimension"] = 2;
        root["Radius"] = m_SGrid.GetDomainRadius();
        { // Description of particles
            YAML::Node node;
            node["Name"] = "particles";
            node["Animated"] = true;
            node["Primitive"] = "Points";
            node["Material"]["Albedo"] = Vector4f(0, 0, 1, 1);
            root["Objects"].push_back(node);
        }
        { // Description of particles
            YAML::Node node;
            node["Name"] = "DEMparticles";
            node["Animated"] = true;
            node["Primitive"] = "Points";
            node["Material"]["Albedo"] = Vector4f(1, 1, 0, 1);
            root["Objects"].push_back(node);
        }
        { // Description of collider
            YAML::Node node;
            node["Name"] = "collider";
            node["Primitive"] = "Points";
            node["Material"]["Albedo"] = Vector4f(.5f, .5f, .5f, 1);
            root["Objects"].push_back(node);
        }
    }

    void Simulation::Export(std::filesystem::path const &dirname, bool initial) const
    {
        { // Export particles
            std::ofstream fout(dirname / "particles.out", std::ios::binary);
            IO::Write(fout, static_cast<std::uint32_t>(m_Particles.size()));
            for (auto const &particle : m_Particles)
            {
                IO::Write(fout, particle.Position.cast<float>().eval());
            }
            for (auto const &particle : m_Particles)
            {
                IO::Write(fout, static_cast<float>(m_SGrid.GetSpacing() / m_SeedingSubFactor / 2.));
            }
        }
        { // Export regular particles
            std::ofstream fout(dirname / "DEMparticles.out", std::ios::binary);
            IO::Write(fout, static_cast<std::uint32_t>(m_DEMParticles.size()));
            for (auto const &p : m_DEMParticles)
            {
                IO::Write(fout, p.Position.cast<float>().eval());
            }
            for (auto const &p : m_DEMParticles)
            {
                IO::Write(fout, static_cast<float>(m_ParticleRadius));
            }
        }
        if (initial)
        { // Export the collider
            std::ofstream fout(dirname / "collider.out", std::ios::binary);
            IO::Write(fout, static_cast<std::uint32_t>(m_ColliderParticles.size()));
            for (auto const &particle : m_ColliderParticles)
            {
                IO::Write(fout, particle.Position.cast<float>().eval());
            }
            for (auto const &particle : m_ColliderParticles)
            {
                IO::Write(fout, static_cast<float>(m_SGrid.GetSpacing() / m_SeedingSubFactor / 2.));
            }
        }
    }

    void Simulation::Save(std::ostream &out) const {}

    void Simulation::Load(std::istream &in) {}

    double Simulation::GetCourantTimeStep() const
    {
        double maxVel = 0;
        for (auto const &particle : m_DEMParticles)
        {
            maxVel = std::max(maxVel, particle.Velocity.norm());
        }
        // printf("%lf\n",maxVel);
        maxVel += std::sqrt(9.8 * m_ParticleRadius * 2);
        return std::min(m_SGrid.GetSpacing() / m_Velocity.GetMaxAbsComponent() * 2, m_ddt * 1000);
    }

    void Simulation::Initialize()
    {
        m_Collider.Finish(m_SGrid);
        CSG::Intersect(m_LevelSet, m_Collider.GetDomainBox());

        SeedParticles();
        ReconstructLevelSet();

        m_NumPartPerCell = m_SeedingSubFactor * m_SeedingSubFactor;
        if (m_DensityCorrectionEnabled)
            m_Pressure.InitRestDensity(m_NumPartPerCell, m_ColliderParticles, m_RestDensity);

        m_ParticleMass = m_DEMDensity * m_ParticleRadius * m_ParticleRadius * std::numbers::pi;
        m_ddt = m_ParticleRadius * std::numbers::pi * std::sqrt(m_DEMDensity / DEMParticle::Young);
        m_CouplingForce.SetZero();
        m_TargetFraction.SetConstant(1.);
        m_DEMForce = DEMForce(m_ParticleRadius);

        std::cerr << "m_single_ratio: " << m_single_ratio << std::endl;
    }

    void Simulation::Advance(double deltaTime)
    {
        MoveDEMParticles(deltaTime);

        AdvectFields(deltaTime);

        CalDensity();

        CalFraction();

        ProjectDensity();

        TransferFromParticlesToGrid();

        ApplyBodyForces(deltaTime);

        ProjectVelocity();

        TransferFromGridToParticles();

        m_LastDeltaTime = deltaTime;
    }

    void Simulation::TransferFromGridToParticles()
    {
        switch (m_Scheme)
        {
        case Scheme::PIC:
            tbb::parallel_for_each(m_Particles.begin(), m_Particles.end(), [&](Particle &particle)
                                   {
                for (int axis = 0; axis < 2; axis++) {
                    particle.Velocity[axis] = BiLerp::Interpolate(m_Velocity[axis], particle.Position);
                } });
            break;
        case Scheme::FLIP:
            tbb::parallel_for_each(m_Particles.begin(), m_Particles.end(), [&](Particle &particle)
                                   {
                for (int axis = 0; axis < 2; axis++) {
                    particle.Velocity[axis] += BiLerp::Interpolate(m_VelDiff[axis], particle.Position);
                    particle.Velocity[axis] = m_BlendingFactor * particle.Velocity[axis]
                        + (1 - m_BlendingFactor) * BiLerp::Interpolate(m_Velocity[axis], particle.Position);
                } });
            break;
        case Scheme::APIC:
            tbb::parallel_for_each(m_Particles.begin(), m_Particles.end(), [&](Particle &particle)
                                   {
                for (int axis = 0; axis < 2; axis++) {
                    particle.Velocity[axis] = BiLerp::Interpolate(m_Velocity[axis], particle.Position);
                    particle.VelocityDrv[axis].setZero();
                    for (auto const [face, weight] : BiLerp::GetGradWtPoints(m_Velocity[axis].GetGrid(), particle.Position)) {
                        particle.VelocityDrv[axis] += m_Velocity[axis].At(face) * weight;
                    }
                } });
            break;
        }
    }

    void Simulation::TransferFromParticlesToGrid()
    {
        m_Velocity.SetZero();
        SGridData<double> weightSum(m_Velocity.GetGrids());

        if (m_Scheme != Scheme::APIC)
        {
            for (auto const &particle : m_Particles)
            {
                for (int axis = 0; axis < 2; axis++)
                {
                    for (auto [face, weight] : BiLerp::GetWtPoints(m_Velocity[axis].GetGrid(), particle.Position))
                    {
                        face = m_Velocity[axis].GetGrid().Clamp(face);
                        m_Velocity[axis][face] += particle.Velocity[axis] * weight;
                        weightSum[axis][face] += weight;
                    }
                }
            }
        }
        else
        {
            for (auto const &particle : m_Particles)
            {
                for (int axis = 0; axis < 2; axis++)
                {
                    for (auto [face, weight] : BiLerp::GetWtPoints(m_Velocity[axis].GetGrid(), particle.Position))
                    {
                        Vector2d const deltaPos = m_Velocity[axis].GetGrid().PositionOf(face) - particle.Position;
                        face = m_Velocity[axis].GetGrid().Clamp(face);
                        m_Velocity[axis][face] += (particle.Velocity[axis] + particle.VelocityDrv[axis].dot(deltaPos)) * weight;
                        weightSum[axis][face] += weight;
                    }
                }
            }
        }

        ParallelForEach(m_Velocity.GetGrids(), [&](int axis, Vector2i const &face)
                        {
            if (weightSum[axis][face]) { m_Velocity[axis][face] /= weightSum[axis][face]; } });

        Extrapolation::Solve(m_Velocity, 0., 3, [&](int axis, Vector2i const &face)
                             { return weightSum[axis][face] != 0; });

        m_VelDiff = m_Velocity;
    }

    void Simulation::AdvectFields(double dt)
    {
        tbb::parallel_for_each(m_Particles.begin(), m_Particles.end(), [&](Particle &particle)
                               { particle.Position = AdvectionSL::Trace<2>(particle.Position, m_Velocity, dt); });

        m_Collider.Enforce(m_Particles);
        ReconstructLevelSet();
    }

    void Simulation::ApplyBodyForces(double dt)
    {
        if (m_GravityEnabled)
        {
            ParallelForEach(m_Velocity[1].GetGrid(), [&](Vector2i const &face)
                            { m_Velocity[1][face] -= 9.8 * dt; });
        }
        double dx = m_SGrid.GetSpacing();
        ParallelForEach(m_Velocity.GetGrids(), [&](int axis, Vector2i const &face)
                        { m_Velocity[axis][face] -= m_CouplingForce[axis][face] * dt / (dx * dx); });
        m_CouplingForce.SetZero();
    }

    void Simulation::ProjectVelocity()
    {
        CalFraction();

        m_Pressure.Project(m_Velocity, m_LevelSet, m_Collider, m_TargetFraction);
        Extrapolation::Solve(m_Velocity, 0., 6, [&](int axis, Vector2i const &face)
                             {
            Vector2i const cell0 = StaggeredGrid::AdjCellOfFace(axis, face, 0);
            Vector2i const cell1 = StaggeredGrid::AdjCellOfFace(axis, face, 1);
            return m_Collider.GetFraction()[axis][face] < 1 && (m_LevelSet[cell0] <= 0 || m_LevelSet[cell1] <= 0); });
        m_Collider.Enforce(m_Velocity);

        // \partial v / \partial t
        ParallelForEach(m_VelDiff.GetGrids(), [&](int axis, Vector2i const &face)
                        { m_VelDiff[axis][face] = m_Velocity[axis][face] - m_VelDiff[axis][face]; });

        // (v /dot /grad) v convection diriv
        // m_ConvectVelocity.SetZero();
        // ParallelForEach(m_ConvectVelocity.GetGrids(), [&](int axis, Vector2i const & face) {
        //     for (int ax = 0; ax < 2; ax++) {
        //         m_ConvectVelocity[axis][face] += (m_Velocity[axis][m_Velocity[axis].GetGrid().Clamp(face + Vector2i::Unit(ax))]
        //                                           - m_Velocity[axis][m_Velocity[axis].GetGrid().Clamp(face -
        //                                           Vector2i::Unit(ax))])
        //             * m_Velocity[ax][face] / m_SGrid.GetSpacing() * 0.5;
        //     }
        // });
    }

    void Simulation::ProjectDensity()
    {
        if (m_DensityCorrectionEnabled)
        {
            ReconstructLevelSet();
            m_Pressure.Correct(m_Particles, m_LevelSet, m_Collider, m_TargetFraction, m_Density);
            m_Collider.Enforce(m_Particles);
            ReconstructLevelSet();
        }
    }

    void Simulation::SeedParticles()
    {
        Extrapolation::Solve(
            m_LevelSet, 1.5 * m_SGrid.GetSpacing(), 1, [&](Vector2i const &cell)
            { return !m_Collider.IsInside(cell); });
        RedistancingFM::Solve(m_LevelSet, 5);

        // A stratified sampling
        auto const rgrid = m_SGrid.CreateRefinedCellGrid(m_SeedingSubFactor);
        double const dx = m_SGrid.GetSpacing();
        double const rad = m_ParticleRadFactor * dx;

        m_Particles.clear();
        ForEach(rgrid, [&](Vector2i const &cell)
                {
            Vector2d const pos = rgrid.PositionOf(cell) + Vector2d::Random() * rgrid.GetSpacing() / 2;
            if (BiLerp::Interpolate(m_Collider.LevelSet, pos) <= 0) {
                m_ColliderParticles.push_back(
                    {
                        .Position = pos,
                    });
            } else if (BiLerp::Interpolate(m_LevelSet, pos) + rad <= 0) {
                m_Particles.push_back(
                    {
                        .Position = pos,
                    });
            } });
    }

    void Simulation::ReconstructLevelSet()
    {
        double const dx = m_SGrid.GetSpacing();
        double const rad = m_ParticleRadFactor * dx;
        m_LevelSet.SetConstant(2 * dx - rad);
        for (auto const &particle : m_Particles)
        {
            for (auto const &cell : BiBSpline<3>::GetPoints(m_LevelSet.GetGrid(), particle.Position))
            {
                if (!m_LevelSet.GetGrid().IsValid(cell))
                    break;
                Vector2d const pos = m_LevelSet.GetGrid().PositionOf(cell);
                m_LevelSet[cell] = std::min(m_LevelSet[cell], (particle.Position - pos).norm() - rad);
            }
        }
        RedistancingFM::Solve(m_LevelSet, 5);

        if (m_SmoothSurfaceEnabled)
        {
            auto const oldLevelSet = m_LevelSet;
            ParallelForEach(m_LevelSet.GetGrid(), [&](Vector2i const &cell)
                            {
                double mean = 0;
                for (int i = 0; i < Grid::GetNumNeighbors(); i++) {
                    Vector2i const nbCell = Grid::NeighborOf(cell, i);
                    mean += oldLevelSet.At(nbCell);
                }
                mean /= Grid::GetNumNeighbors();
                if (mean < m_LevelSet[cell]) m_LevelSet[cell] = mean; });
            RedistancingFM::Solve(m_LevelSet, 5);
        }
    }

    void Simulation::CacheNeighborHoods()
    {
        ParallelForEach(m_DEMGrid.GetGrid(), [&](Vector2i const &cell)
                        { m_DEMGrid[cell].clear(); });
        for (const auto &p : m_DEMParticles)
        {
            Vector2i lower = m_DEMGrid.GetGrid().Clamp(
                m_DEMGrid.GetGrid().CalcLower<1>(p.Position + Vector2d::Ones() * m_SGrid.GetSpacing() * 0.5));
            m_DEMGrid[lower].push_back(p);
        }
    }

    void Simulation::MoveDEMParticles(double dt)
    {
        double deltaTime = dt;

        while (1)
        {
            double maxVel = tbb::parallel_reduce(
                tbb::blocked_range<size_t>(0, m_DEMParticles.size()),
                0.0,
                [&](const tbb::blocked_range<size_t> &range, double local_max) -> double
                {
                    for (size_t i = range.begin(); i != range.end(); ++i)
                    {
                        local_max = std::max(local_max, m_DEMParticles[i].Velocity.norm());
                    }
                    return local_max;
                },
                [](double a, double b) -> double
                { return std::max(a, b); });
            // printf("%lf\n",maxVel);
            maxVel += std::sqrt(9.8 * m_ParticleRadius * 2);
            double ddt = std::min(m_ParticleRadius * 2. / maxVel, m_ddt);

            if (dt < ddt)
                break;
            else
                dt -= ddt;

            MoveDEMParticlesSplit(ddt, deltaTime);
        }
        MoveDEMParticlesSplit(dt, deltaTime);
    }

    void Simulation::MoveDEMParticlesSplit(double ddt, double dt)
    {
        tbb::parallel_for_each(m_DEMParticles.begin(), m_DEMParticles.end(), [&](DEMParticle &particle)
                               {
                                   double den = BiLerp::Interpolate(m_Density, particle.Position);
                                   particle.SaturateRate = std::clamp(den / (1 - den), 0., 1.); //
                               });
        CacheNeighborHoods();
        CalCoupling(ddt, dt);
        ApplyDEMForces(ddt, dt);
        TransferCouplingForces(ddt, dt);
        m_Collider.Enforce(m_DEMParticles, ddt, m_ParticleRadius);
        tbb::parallel_for_each(m_DEMParticles.begin(), m_DEMParticles.end(), [&](DEMParticle &particle)
                               {
            particle.Velocity += particle.AccVelocity * ddt;
            particle.Position += particle.Velocity * ddt; });
    }

    void Simulation::ApplyDEMForces(double ddt, double dt)
    {
        // Note: Do not apply external forces to boundary particles
        tbb::parallel_for_each(m_DEMParticles.begin(), m_DEMParticles.end(), [&](DEMParticle &particle)
                               {
            auto ac = particle.AccVelocity;
            particle.AccVelocity = Vector2d::Zero();
            particle.AccVelocity[1] -= 9.8;
            if (BiLerp::Interpolate(m_LevelSet, particle.Position) <= 0)
            {
                particle.AccVelocity += (particle.CouplingForce + m_DEMForce.getForceSum(m_DEMGrid, particle) + (BiLerp::Interpolate(m_VelDiff, particle.Position) / m_LastDeltaTime
                                                                                                                 // + BiLerp::Interpolate(m_ConvectVelocity, particle.Position)
                                                                                                                 ) * m_ParticleVolume *
                                                                                                                    0.5) /
                                        (m_ParticleMass + m_ParticleVolume * 0.5);
            }
            else 
                particle.AccVelocity += (particle.CouplingForce + m_DEMForce.getForceSum(m_DEMGrid, particle)) / (m_ParticleMass); });
    }

    void Simulation::CalCoupling(double ddt, double dt)
    {
        tbb::parallel_for_each(m_DEMParticles.begin(), m_DEMParticles.end(), [&](DEMParticle &particle)
                               {
            if (BiLerp::Interpolate(m_LevelSet, particle.Position) <= 0) {
                particle.CouplingForce -= ((BiLerp::Interpolate(m_Pressure.GetGradPressure1(), particle.Position)
                                            + BiLerp::Interpolate(m_Pressure.GetGradPressure2(), particle.Position)
                                        )
                                           / m_LastDeltaTime)
                    * m_ParticleVolume;
                particle.CouplingForce += (BiLerp::Interpolate(m_Velocity, particle.Position) - particle.Velocity)
                    * (BiLerp::Interpolate(m_Velocity, particle.Position) - particle.Velocity).norm() * 3 * std::numbers::pi
                    * m_ParticleRadius * m_ViscosityCoeff;
                particle.CouplingForce[1] += m_ParticleVolume*9.8;
            } });
    }
    void Simulation::TransferCouplingForces(double ddt, double dt)
    {
        SGridData<double> weightSum(m_CouplingForce.GetGrids());
        SGridData<double> m_RestCouplingForce(m_CouplingForce.GetGrids());
        for (auto const &particle : m_DEMParticles)
        {
            for (int axis = 0; axis < 2; axis++)
            {
                for (auto [face, weight] : BiLerp::GetWtPoints(m_RestCouplingForce[axis].GetGrid(), particle.Position))
                {
                    face = m_RestCouplingForce[axis].GetGrid().Clamp(face);
                    m_RestCouplingForce[axis][face] += particle.CouplingForce[axis] * weight;
                    weightSum[axis][face] += weight;
                }
            }
        }

        ParallelForEach(m_CouplingForce.GetGrids(), [&](int axis, Vector2i const &face)
                        {
            if (weightSum[axis][face]) { m_CouplingForce[axis][face] += m_RestCouplingForce[axis][face] * ddt / (dt); } });

        tbb::parallel_for_each(m_DEMParticles.begin(), m_DEMParticles.end(), [&](DEMParticle &particle)
                               { particle.CouplingForce = Vector2d::Zero(); });
    }

    void Simulation::CalFraction()
    {
        m_TargetFraction.SetConstant(1.);
        CacheNeighborHoods();
        double frac = m_ParticleVolume / m_SGrid.GetSpacing() / m_SGrid.GetSpacing();
        for (const auto &particle : m_DEMParticles)
        {
            for (auto [face, weight] : BiLerp::GetWtPoints(m_TargetFraction.GetGrid(), particle.Position))
            {
                face = m_TargetFraction.GetGrid().Clamp(face);
                m_TargetFraction[face] -= weight * frac;
            }
        }

        ParallelForEach(m_TargetFraction.GetGrid(), [&](const Vector2i &cell)
                        { m_TargetFraction[cell] = std::max(m_TargetFraction[cell], 1 - 0.907); });
    }

    void Simulation::CalDensity()
    {
        m_Density = m_RestDensity;
        for (auto const &particle : m_Particles)
        {
            for (auto const [cell, weight] : BiLerp::GetWtPoints(m_Density.GetGrid(), particle.Position))
            {
                m_Density[m_Density.GetGrid().Clamp(cell)] += weight / m_NumPartPerCell;
            }
        }
    }
} // namespace Pivot
