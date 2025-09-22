#pragma once

#include "Collider.h"
#include "DEMForce.h"
#include "Pressure.h"

namespace Pivot
{
    class Simulation
    {
    private:
        friend class SimBuilder;

    public:
        enum class Scheme
        {
            PIC,
            FLIP,
            APIC
        };
        enum class Scene
        {
            Falling,
            BigBall
        };

    public:
        explicit Simulation(StaggeredGrid const &sgrid, double particleRadius);

        void Describe(YAML::Node &root) const;
        void Export(std::filesystem::path const &dirname, bool initial = false) const;
        void Save(std::ostream &out) const;
        void Load(std::istream &in);

        void Initialize();
        void Advance(double deltaTime);

        void TransferFromGridToParticles();
        void TransferFromParticlesToGrid();
        void AdvectFields(double dt);
        void ApplyBodyForces(double dt);
        void ProjectVelocity();
        void ProjectDensity();

        void SeedParticles();
        void ReconstructLevelSet();

        void SetTime(double time) { m_Time = time; }
        auto GetTime() const { return m_Time; }

        double GetCourantTimeStep() const;

        void CacheNeighborHoods();
        void CalFraction();
        void CalDensity();

        void MoveDEMParticles(double dt);
        void MoveDEMParticlesSplit(double ddt, double dt);
        void ApplyDEMForces(double ddt, double dt);
        void CalCoupling(double ddt, double dt);
        void TransferCouplingForces(double ddt, double dt);

    private:
        // Parameters for basic simulation
        double m_Time;
        Scene m_Scene;
        double delay_water = 0.0;
        // Grid-based data structures
        StaggeredGrid m_SGrid;
        Collider m_Collider;
        Pressure m_Pressure;
        SGridData<double> m_Velocity;
        SGridData<double> m_ConvectVelocity;
        GridData<double> m_LevelSet;
        SGridData<double> m_VelDiff; // Used for FLIP
        GridData<double> m_RestDensity;
        GridData<double> m_Density;
        // Particle-based data structures
        std::vector<Particle> m_ColliderParticles;
        std::vector<Particle> m_Particles;

        // Parameters for the scheme
        Scheme m_Scheme = Scheme::PIC;
        double m_BlendingFactor = 0.95; // Used for FLIP
        double m_LastDeltaTime = 1e6;   // very large
        double m_ViscosityCoeff = 0.44;
        // Parameters for particles
        int m_SeedingSubFactor = 3;
        int m_NumPartPerCell;
        double m_ParticleRadFactor = 1.01 * std::numbers::sqrt2 / 2;
        // Boolean configurations
        bool m_DensityCorrectionEnabled = true;
        bool m_SmoothSurfaceEnabled = true;
        bool m_GravityEnabled = true;

        // DEM structure

        double m_ParticleRadius;
        double m_SupportRadius;
        double m_ParticleMass;
        double m_ParticleVolume;

        double m_single_ratio;

        DEMForce m_DEMForce;

        GridData<std::vector<DEMParticle>> m_DEMGrid;

        std::vector<DEMParticle> m_DEMParticles;

        double m_DEMDensity = 2.5;
        double m_ddt;

        // coupling structure
        SGridData<double> m_CouplingForce;
        GridData<double> m_TargetFraction;
    };
} // namespace Pivot
