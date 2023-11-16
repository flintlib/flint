/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fexpr.h"
#include "fexpr_builtin.h"

static void
traverse(fexpr_vec_t nodes, const fexpr_t expr)
{
    slong i, nargs;
    fexpr_t view;

    if (fexpr_is_integer(expr))
        return;

    if (fexpr_is_arithmetic_operation(expr))
    {
        nargs = fexpr_nargs(expr);

        fexpr_view_arg(view, expr, 0);
        for (i = 0; i < nargs; i++)
        {
            traverse(nodes, view);
            fexpr_view_next(view);
        }
        return;
    }

    if (fexpr_is_builtin_call(expr, FEXPR_Pow) && (fexpr_nargs(expr) == 2))
    {
        fexpr_t base, exp;

        fexpr_view_arg(base, expr, 0);
        fexpr_view_arg(exp, expr, 1);

        if (fexpr_is_integer(exp))
        {
            traverse(nodes, base);
            return;
        }
    }

    fexpr_vec_insert_unique(nodes, expr);
}

void
fexpr_arithmetic_nodes(fexpr_vec_t nodes, const fexpr_t expr)
{
    fexpr_vec_set_length(nodes, 0);
    traverse(nodes, expr);
}
