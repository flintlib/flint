/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

truth_t
gr_mat_is_hessenberg(const gr_mat_t mat, gr_ctx_t ctx)
{
    slong i, r, c;
    truth_t is_zero, row_is_zero;

    r = gr_mat_nrows(mat, ctx);
    c = gr_mat_ncols(mat, ctx);

    if (r != c)
        return T_FALSE;

    is_zero = T_TRUE;

    for (i = 2; i < c && is_zero != T_FALSE; i++)
    {
        row_is_zero = _gr_vec_is_zero(mat->rows[i], i - 1, ctx);
        is_zero = truth_and(is_zero, row_is_zero);
    }

    return is_zero;
}
