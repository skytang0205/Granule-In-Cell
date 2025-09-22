#pragma once

#include "Common.h"

namespace Pivot {
    class Particle {
    public:
        Vector2d Position;
        Vector2d Velocity = Vector2d::Zero();

        std::array<Vector2d, 2> VelocityDrv = { Vector2d::Zero().eval(), Vector2d::Zero().eval() }; // Used for APIC
    };

    class DEMParticle : public Particle {
    public:
        Vector2d AccVelocity   = Vector2d::Zero();
        Vector2d CouplingForce = Vector2d::Zero();

        double SaturateRate = 0.;

        static constexpr double Young     = 1e6;
        static constexpr double Poisson   = .3;
        static constexpr double FricAngle = .5;
    };
} // namespace Pivot
