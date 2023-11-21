/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

int
fexpr_contains(const fexpr_t expr, const fexpr_t x)
{
    fexpr_t func, arg;
    slong expr_size, x_size, i, nargs;

    expr_size = fexpr_size(expr);
    x_size = fexpr_size(x);

    if (expr_size < x_size)
        return 0;

    if (expr_size == x_size)
        return _mpn_equal(expr->data, x->data, x_size);

    nargs = fexpr_nargs(expr);

    if (nargs < 0)
        return 0;

    fexpr_view_func(func, expr);
    if (fexpr_contains(func, x))
        return 1;

    if (nargs <= 0)
        return 0;

    fexpr_view_arg(arg, expr, 0);
    for (i = 0; i < nargs; i++)
    {
        if (fexpr_contains(arg, x))
            return 1;

        if (i < nargs - 1)
            fexpr_view_next(arg);
    }

    return 0;
}
