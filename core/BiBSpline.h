#pragma once

#include "GridData.h"

namespace Pivot {
    template<int Order> class BiBSpline {
    public:
        static std::array<Vector2i, (Order + 1) * (Order + 1)> GetPoints(Grid const & grid, Vector2d const & pos) {
            std::array<Vector2i, (Order + 1) * (Order + 1)> pts;
            Vector2i const                                  lower = grid.CalcLower<Order>(pos);
            for (int i = 0; i <= Order; i++) {
                for (int j = 0; j <= Order; j++) { pts[i * Order + i + j] = lower + Vector2i(i, j); }
            }
            return pts;
        }
    };
}