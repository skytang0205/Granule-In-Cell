#include "VolumeSampler.h"

#include "PoissonDiskSampler.h"

namespace Pivot {
    std::vector<Vector2d> VolumeSampler::Sample(Surface const & surface, bool if_Poission) {
        std::vector<Vector2d> samples;

        if (if_Poission) {
            PoissonDiskSampler sampler(radius, surface);
            samples = sampler.generatePoints(surface);
        } else {
            auto const [minCorner, maxCorner] = surface.GetCornersOfAABB();
            auto const minCell                = CellOf(minCorner);
            auto const maxCell                = CellOf(maxCorner);
            for (auto const & cell : IterateOverCells(minCell, maxCell)) {
                Vector2d const center = CenterOf(cell);
                for (int axis = 0; axis < 2; ++axis) {
                    Vector2d const faceCenter0 = center - Vector2d::Unit(axis) * m_CellSize[axis] * .5;
                    if (surface.Surrounds(faceCenter0)) samples.push_back(faceCenter0);
                    if (cell[axis] == maxCell[axis]) {
                        Vector2d const faceCenter1 = center + Vector2d::Unit(axis) * m_CellSize[axis] * .5;
                        if (surface.Surrounds(faceCenter1)) samples.push_back(faceCenter1);
                    }
                }
            }
        }
        std::cout << samples.size() << std::endl;
        return samples;
    }

    std::vector<Vector2i> VolumeSampler::IterateOverCells(Vector2i const & minCell, Vector2i const & maxCell) const {
        std::vector<Vector2i> cells;
        cells.reserve((maxCell.x() - minCell.x() + 1) * (maxCell.y() - minCell.y() + 1));
        for (int x = minCell.x(); x <= maxCell.x(); ++x) {
            for (int y = minCell.y(); y <= maxCell.y(); ++y) { cells.emplace_back(x, y); }
        }
        return cells;
    }
} // namespace Pivot