#include "SimBuilder.h"
#include "VolumeSampler.h"

#include "CSG.h"

namespace Pivot
{
    std::unique_ptr<Simulation> SimBuilder::Build(SimBuildOptions const &options)
    {
        std::unique_ptr<Simulation> simulation;
        switch (options.Scene)
        {
        case Simulation::Scene::Falling:
            simulation = BuildFalling(options);
            break;
        case Simulation::Scene::BigBall:
            simulation = BuildBigBall(options);
            break;
        }
        simulation->m_Scene = options.Scene;
        return simulation;
    }

    std::unique_ptr<Simulation> SimBuilder::BuildFalling(SimBuildOptions const &options)
    {
        constexpr double length = 1.;
        constexpr int bw = 2;
        int const scale = options.Scale < 0 ? 128 : options.Scale;
        double const radius = options.ParticleRadius < 0 ? .5 / double(scale) : options.ParticleRadius / double(scale);
        StaggeredGrid sgrid(2, length / (scale - bw * 2), Vector2i(1, 2) * scale);
        auto sim = std::make_unique<Simulation>(sgrid, radius);
        AddDEMParticles(sim.get(), ImplicitSphere(Vector2d::Unit(1) * length * -0.2, length * .15), true);
        CSG::Union(sim->m_LevelSet, ImplicitBox(-Vector2d(.5, 1.05) * length, Vector2d(1, 0.5) * length));
        return sim;
    }

    std::unique_ptr<Simulation> SimBuilder::BuildBigBall(SimBuildOptions const &options)
    {
        constexpr double length = 1.;
        constexpr int bw = 2;
        int const scale = options.Scale < 0 ? 128 : options.Scale;
        double const radius = options.ParticleRadius < 0 ? .5 / double(scale) : options.ParticleRadius / double(scale);
        StaggeredGrid sgrid(2, length / (scale - bw * 2), Vector2i(1, 2) * scale);
        auto sim = std::make_unique<Simulation>(sgrid, radius);
        CSG::Union(sim->m_LevelSet, ImplicitSphere(Vector2d::Unit(1) * length * -0.3, length * .35));
        AddDEMParticles(sim.get(), ImplicitSphere(Vector2d::Unit(1) * length * -0.3, length * .35), false);
        return sim;
    }

    void SimBuilder::AddDEMParticles(
        Simulation *simulation, Surface const &surface, bool if_Poission, std::function<Vector2d(Vector2d const &)> velocity)
    {
        VolumeSampler sampler(simulation->m_ParticleRadius * 2.);
        auto const positions = sampler.Sample(surface, if_Poission);
        simulation->m_DEMParticles.reserve(simulation->m_DEMParticles.size() + positions.size());
        for (auto const &pos : positions)
        {
            simulation->m_DEMParticles.push_back(
                {{
                    .Position = pos,
                    .Velocity = velocity ? velocity(pos) : Vector2d::Zero(),
                }});
        }
    }
} // namespace Pivot
