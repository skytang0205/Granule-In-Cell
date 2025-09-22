#pragma once

#include "Surface.h"

namespace Pivot {
    class VolumeSampler {
    public:
        VolumeSampler(double latticeConstant):
            radius(latticeConstant), m_CellSize { latticeConstant, latticeConstant * std::numbers::sqrt3 } {}

        std::vector<Vector2d> Sample(Surface const & surface, bool if_Poission = false);

    private:
        Vector2i CellOf(Vector2d const & pos) const {
            return (pos.array() / m_CellSize - c_Offset).floor().template cast<int>().matrix();
        }
        Vector2d CenterOf(Vector2i const & cell) const {
            return (cell.template cast<double>() + Vector2d::Constant(.5 + c_Offset)).cwiseProduct(m_CellSize.matrix());
        }

        std::vector<Vector2i> IterateOverCells(Vector2i const & minIndex, Vector2i const & maxIndex) const;

    private:
        static constexpr double c_Offset = .02407112028; // A magic number to avoid conflicts

        double const radius;

        Array2d const m_CellSize;
    };
}