/*
    Copyright (C) 2024 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

truth_t
gr_mat_contains_zero_mat(const gr_mat_t mat, gr_ctx_t ctx)
{
    truth_t eq, this_eq;
    slong i, r, c;

    r = gr_mat_nrows(mat, ctx);
    c = gr_mat_ncols(mat, ctx);

    if (r == 0 || c == 0)
        return T_TRUE;

    eq = T_TRUE;

    for (i = 0; i < r; i++)
    {
        this_eq = _gr_vec_contains_zero_vec(mat->rows[i], c, ctx);

        if (this_eq == T_FALSE)
            return T_FALSE;
    }

    return eq;
}
