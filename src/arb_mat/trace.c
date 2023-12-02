/*
    Copyright (C) 2015 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_trace(arb_t trace, const arb_mat_t mat, slong prec)
{
    slong i;

    if (!arb_mat_is_square(mat))
    {
        flint_throw(FLINT_ERROR, "arb_mat_trace: a square matrix is required!\n");
    }

    if (arb_mat_is_empty(mat))
    {
        arb_zero(trace);
        return;
    }

    arb_set(trace, arb_mat_entry(mat, 0, 0));
    for (i = 1; i < arb_mat_nrows(mat); i++)
    {
        arb_add(trace, trace, arb_mat_entry(mat, i, i), prec);
    }
}
