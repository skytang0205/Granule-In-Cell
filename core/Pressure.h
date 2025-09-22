#pragma once

#include "Collider.h"

namespace Pivot {
    class Pressure {
    public:
        explicit Pressure(StaggeredGrid const & sgrid);

    public:
        void Project(
            SGridData<double> &      velocity,
            GridData<double> const & levelSet,
            Collider const &         collider,
            GridData<double> const & fraction);
        void Correct(
            std::vector<Particle> &  particles,
            GridData<double> const & levelSet,
            Collider const &         collider,
            GridData<double> const & fraction,
            GridData<double> const & density);

        void InitRestDensity(int numPartPerCell, std::vector<Particle> const & particles, GridData<double> & density);

        GridData<double> const &  GetPressure() const { return m_Pressure; }
        SGridData<double> const & GetGradPressure1() const { return m_GradPressure1; }
        SGridData<double> const & GetGradPressure2() const { return m_GradPressure2; }

    private:
        void BuildProjectionMatrix(
            SGridData<double> const & velocity,
            GridData<double> const &  levelSet,
            Collider const &          collider,
            GridData<double> const &  fraction);
        void
            BuildCorrectionMatrix(GridData<double> const & density, Collider const & collider, GridData<double> const & fraction);

        void SetUnKnowns(GridData<double> const & levelSet);

        void SolveLinearSystem();

        void ApplyProjection(SGridData<double> & velocity, GridData<double> const & levelSet, Collider const & collider);
        void ApplyCorrection(std::vector<Particle> & particles, Collider const & collider);

    private:
        GridData<int>    m_Grid2Mat;
        std::vector<int> m_Mat2Grid;

        SparseMatrix<double, RowMajor> m_MatL; // A matrix of the Laplacian operator

        VectorXd m_RdP; // reduced pressure
        VectorXd m_Rhs;

        int               m_NumPartPerCell;
        GridData<double>  m_Pressure;
        SGridData<double> m_GradPressure1;
        SGridData<double> m_GradPressure2;
        SGridData<double> m_DeltaPos;
    };
} // namespace Pivot
