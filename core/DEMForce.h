#pragma once

#include "Common.h"
#include "Particle.h"

namespace Pivot
{
    class QuadraticBezierCoeff
    {
    private:
        double a, b, c, d, e;
        double py0, py1, px1, py2;

    public:
        QuadraticBezierCoeff(
            const double py0 = 0., const double py1 = 0., const double px1 = 0., const double py2 = 0.) : py0(py0), py1(py1), px1(px1), py2(py2)
        {
            a = (px1 + 1.f) / 4.f;
            b = -2.f * a;
            c = b * b;
            d = -4.f * (1.f + b - px1);
            e = 2.f * (1.f + b - px1);
        }

        double rx2t(const double sr) const { return (b + std::sqrt(c + d * (px1 - sr))) / e; }

        double calculate(const double sr) const
        {
            if (sr < 0.f)
                return py0;

            if (sr >= 1.f)
                return py2;

            if (sr <= px1)
            {
                const double t = sr / px1;
                const double omt = 1. - t;
                return omt * omt * py0 + 2 * t * omt * py1 + t * t * py1;
            }
            else
            {
                const double t = rx2t(sr);
                const double omt = 1. - t;
                return omt * omt * py1 + 2 * t * omt * py1 + t * t * py2;
            }
        }
    };

    class DEMForce
    {
    public:
        DEMForce(double radius) : _radius(radius)
        {
            volume_liquid_bridge = 4. / 3. * std::numbers::pi * std::pow(radius, 3.) * 0.01 * 0.01;
            d_rupture = std::pow(volume_liquid_bridge, 1. / 3.) + 0.1 * std::pow(volume_liquid_bridge, 2. / 3.);
            c0 = 0.0;
            cmc = 1.;
            cmcp = 0.1;
            csat = -1.5;
            G = QuadraticBezierCoeff(c0, cmc, cmcp, csat);
            surface_tensor_cof = 0.007;
            printf("drupture: %lf\n", d_rupture / radius);
        }
        double _tan_fricangle() const { return std::tan(DEMParticle::FricAngle); }

        double _K_norm() const { return DEMParticle::Young * _radius; }

        double _K_tang() const { return _K_norm() * DEMParticle::Poisson; }

        Vector2d ComputeDemCapillaryForces(Vector2d dij, double sr) const
        {
            Vector2d f = Vector2d::Zero();
            double sij = dij.norm() - 2 * _radius;
            if (sij < d_rupture && sij > 0.000001)
            {
                Vector2d n = dij.normalized();

                double coeff_c = G.calculate(sr) * surface_tensor_cof;

                double d = -sij + std::sqrt(sij * sij + volume_liquid_bridge / (std::numbers::pi * _radius));

                f = -2. * std::numbers::pi * coeff_c * _radius / (1. + sij / (2. * d)) * -n;
            }
            return f;
        }

        Vector2d getForce(const DEMParticle &pi, const DEMParticle &pj) const
        {
            Vector2d dij = pj.Position - pi.Position;
            Vector2d vij = pj.Velocity - pi.Velocity;
            double sr = (pi.SaturateRate + pj.SaturateRate) * .5;

            Vector2d force = Vector2d::Zero();

            {
                // Compute DEM Force
                Vector2d f = Vector2d::Zero();
                double dist = dij.norm();
                double penetration_depth = 2 * _radius - dist;
                if (penetration_depth > 0. && dist > 0.0001 * _radius)
                {
                    Vector2d n = dij.normalized();
                    double dot_epslion = vij.dot(n);
                    Vector2d vij_tangential = vij - dot_epslion * n;

                    Vector2d normal_force = _K_norm() * penetration_depth * n;
                    Vector2d shear_force = -_K_tang() * vij_tangential;

                    double max_fs = normal_force.norm() * _tan_fricangle();
                    double shear_force_norm = shear_force.norm();

                    if (shear_force_norm > max_fs)
                    {
                        shear_force = shear_force * max_fs / shear_force_norm;
                    }
                    f = -normal_force - shear_force;
                }
                force += f;
            }
            return force + ComputeDemCapillaryForces(dij, sr);
        }

        Vector2d getForceSum(GridData<std::vector<DEMParticle>> const &grid, const DEMParticle &particle) const
        {
            Vector2d f = Vector2d::Zero();
            Vector2i lower = grid.GetGrid().Clamp(grid.GetGrid().CalcLower<1>(particle.Position));
            int range = int(_radius * grid.GetGrid().GetInvSpacing() + 2.) * 2;
            Vector2i size = grid.GetGrid().GetSize();
            for (int i = std::max(0, lower.x() - range); i < std::min(size.x(), lower.x() + range); i++)
            {
                for (int j = std::max(0, lower.y() - range); j < std::min(size.y(), lower.y() + range); j++)
                {
                    for (auto const &nb : grid[Vector2i(i, j)])
                    {
                        if ((nb.Position - particle.Position).squaredNorm() < 6 * _radius * _radius)
                            f += getForce(particle, nb);
                    }
                }
            }
            return f;
        }

        Vector2d getForceSum(std::vector<DEMParticle> const &m_Particles, const DEMParticle &p0) const
        {
            Vector2d f = Vector2d::Zero();
            for (auto const &p1 : m_Particles)
            {
                if (((p0.Position - p1.Position).squaredNorm() < 4 * _radius * _radius))
                {
                    f = f + getForce(p0, p1);
                }
            }
            return f;
        }

    private:
        double _radius;
        double volume_liquid_bridge, d_rupture;
        double c0, cmc, cmcp, csat, sr, surface_tensor_cof;
        QuadraticBezierCoeff G;
    };
} // namespace Pivot