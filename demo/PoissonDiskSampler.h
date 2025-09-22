#pragma once

#include "Common.h"
#include "ImplicitSurface.h"
#include "Surface.h"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <vector>

namespace Pivot {
    class PoissonDiskSampler {
    public:
        PoissonDiskSampler(double radius, Surface const & surface): radius(radius) {
            auto const [minCorner, maxCorner] = surface.GetCornersOfAABB();
            low_bound                         = minCorner;
            up_bound                          = maxCorner;
            gridCellSize                      = radius / std::numbers::sqrt2;
            grid_size                         = ((up_bound - low_bound) / gridCellSize).cast<int>();
            int size                          = 1;
            for (int i = 0; i < 2; i++) size *= grid_size(i);
            grid.resize(size, -1);
        }

        Vector2d generateNewPoint(Vector2d point) {
            static std::random_device               rd;
            static std::mt19937                     gen(rd());
            static std::uniform_real_distribution<> dis(0.0, 1.0);

            double   angle    = dis(gen) * 2 * std::numbers::pi;
            double   dist     = radius + dis(gen) * radius;
            Vector2d newPoint = { std::cos(angle), std::sin(angle) };
            newPoint          = newPoint * dist + point;
            return newPoint;
        }

        std::vector<Vector2d> generatePoints(Surface const & surface) {
            std::vector<Vector2d> samples;
            std::vector<Vector2d> processList;
            Vector2d              initialPoint = (up_bound + low_bound) * 0.5;
            processList.push_back(initialPoint);
            addToGrid(initialPoint);

            while (! processList.empty()) {
                int      idx   = rand() % processList.size();
                Vector2d point = processList[idx];
                bool     found = false;

                for (int i = 0; i < 300; ++i) {
                    Vector2d newPoint = generateNewPoint(point);
                    if (isInBounds(newPoint) && surface.Surrounds(newPoint) && isFarEnough(newPoint)) {
                        // samples.push_back(newPoint);
                        processList.push_back(newPoint);
                        addToGrid(newPoint);
                        found = true;
                        break;
                    }
                    // if(!isInBounds(newPoint)) printf("not in bounds\n");
                    // if (!isFarEnough(newPoint)) printf("not far\n");
                }
                if (! found) { processList.erase(processList.begin() + idx); }
                printf("%zd %zd %zd\n", grid.size(), samples.size(), processList.size());
                if (gridPoints.size() > grid.size() * 10) break;
            }

            for (const auto & idx : grid) {
                if (idx != -1) { samples.push_back(gridPoints[idx]); }
            }
            return samples;
        }

        bool isInBounds(const Vector2d & point) const {
            Vector2i grid_idx = ((point - low_bound) / gridCellSize).cast<int>();
            return grid_idx.x() >= 0 && grid_idx.x() < grid_size.x() && grid_idx.y() >= 0 && grid_idx.y() < grid_size.y();
        }

    private:
        double                radius, gridCellSize;
        Vector2d              boundary, center, low_bound, up_bound;
        Vector2i              grid_size;
        std::vector<int>      grid;
        std::vector<Vector2d> gridPoints;

        bool isFarEnough(const Vector2d & point) const {
            Vector2i grid_idx = ((point - low_bound) / gridCellSize).cast<int>();

            for (int y = std::max(0, grid_idx.y() - 3); y <= std::min(grid_size.y() - 1, grid_idx.y() + 3); ++y) {
                for (int x = std::max(0, grid_idx.x() - 3); x <= std::min(grid_size.x() - 1, grid_idx.x() + 3); ++x) {
                    int i   = y * grid_size.x() + x;
                    int idx = grid[i];
                    if (idx != -1) {
                        Vector2d p = gridPoints[idx];
                        if ((p - point).norm() < radius) { return false; }
                    }
                }
            }
            return true;
        }

        void addToGrid(const Vector2d & point) {
            Vector2i grid_idx                                 = ((point - low_bound) / gridCellSize).cast<int>();
            grid[grid_idx.y() * grid_size.x() + grid_idx.x()] = gridPoints.size();
            gridPoints.push_back(point);
        }
    };
} // namespace Pivot