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
fexpr_num_leaves(const fexpr_t expr)
{
    if (fexpr_is_atom(expr))
    {
        return 1;
    }
    else
    {
        fexpr_t func, arg;
        slong i, leaves, nargs;

        fexpr_view_func(func, expr);
        leaves = fexpr_num_leaves(func);

        nargs = fexpr_nargs(expr);
        *arg = *func;

        for (i = 0; i < nargs; i++)
        {
            fexpr_view_next(arg);
            leaves += fexpr_num_leaves(arg);
        }

        return leaves;
    }
}
