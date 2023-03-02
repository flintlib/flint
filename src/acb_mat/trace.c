/*
    Copyright (C) 2015 Arb authors

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_trace(acb_t trace, const acb_mat_t mat, slong prec)
{
    slong i;

    if (!acb_mat_is_square(mat))
    {
        flint_printf("acb_mat_trace: a square matrix is required!\n");
        flint_abort();
    }

    if (acb_mat_is_empty(mat))
    {
        acb_zero(trace);
        return;
    }

    acb_set(trace, acb_mat_entry(mat, 0, 0));
    for (i = 1; i < acb_mat_nrows(mat); i++)
    {
        acb_add(trace, trace, acb_mat_entry(mat, i, i), prec);
    }
}
