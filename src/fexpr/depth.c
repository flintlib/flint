/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

slong
fexpr_depth(const fexpr_t expr)
{
    if (fexpr_is_atom(expr))
    {
        return 1;
    }
    else
    {
        fexpr_t func, arg;
        slong i, depth, d, nargs;

        fexpr_view_func(func, expr);
        depth = fexpr_depth(func);

        nargs = fexpr_nargs(expr);
        *arg = *func;

        for (i = 0; i < nargs; i++)
        {
            fexpr_view_next(arg);
            d = fexpr_depth(arg);
            depth = FLINT_MAX(depth, d);
        }

        return depth + 1;
    }
}
