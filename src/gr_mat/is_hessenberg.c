/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

truth_t
gr_mat_is_hessenberg(const gr_mat_t mat, gr_ctx_t ctx)
{
    gr_method_vec_predicate is_zero = GR_VEC_PREDICATE(ctx, VEC_IS_ZERO);
    slong ar, ac, i, sz = ctx->sizeof_elem;
    truth_t eq, this_eq;

    ar = gr_mat_nrows(mat, ctx);
    ac = gr_mat_ncols(mat, ctx);

    if (ar <= 2 || ac == 0)
        return T_TRUE;

    eq = T_TRUE;

    for (i = 2; i < ar; i++)
    {
        this_eq = is_zero(GR_MAT_ENTRY(mat, i, 0, sz), FLINT_MIN(i - 1, ac), ctx);

        if (this_eq == T_FALSE)
            return T_FALSE;

        if (this_eq == T_UNKNOWN)
            eq = T_UNKNOWN;
    }

    return eq;
}
