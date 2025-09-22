#pragma once

#include "GridData.h"

namespace Pivot {
    class FiniteDiff {
    public:
        template<typename Type> static Type CalcFirstDrv(GridData<Type> const & grData, Vector2i const & coord, int axis) {
            double const invDx = grData.GetGrid().GetInvSpacing();
            auto const   f     = [&](int i) -> Type { return grData[coord + Vector2i::Unit(axis) * i]; };
            if (coord[axis] == 0) {
                return CalcFirstForward(f(0), f(+1), f(+2)) * invDx * .5;
            } else if (coord[axis] + 1 == grData.GetGrid().GetSize()[axis]) {
                return CalcFirstBackward(f(-2), f(-1), f(0)) * invDx * .5;
            } else {
                return CalcFirstCentral(f(-1), f(+1)) * invDx * .5;
            }
        }

    private:
        template<typename Type> static Type CalcFirstCentral(Type const & f0, Type const & f2) { return f2 - f0; }
        template<typename Type> static Type CalcFirstForward(Type const & f0, Type const & f1, Type const & f2) {
            return 4 * f1 - 3 * f0 - f2;
        }
        template<typename Type> static Type CalcFirstBackward(Type const & f0, Type const & f1, Type const & f2) {
            return 3 * f2 - 4 * f1 + f0;
        }
    };
} // namespace Pivot
