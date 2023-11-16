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
gr_mat_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx)
{
    gr_method_vec_predicate is_zero = GR_VEC_PREDICATE(ctx, VEC_IS_ZERO);
    gr_method_unary_predicate is_neg_one = GR_UNARY_PREDICATE(ctx, IS_NEG_ONE);
    truth_t eq, this_eq;
    slong i, ar, ac, sz;

    ar = gr_mat_nrows(mat, ctx);
    ac = gr_mat_ncols(mat, ctx);

    if (ar == 0 || ac == 0)
        return T_TRUE;

    sz = ctx->sizeof_elem;

    eq = T_TRUE;

    for (i = 0; i < ar; i++)
    {
        /* upper triangular */
        if (i >= 1)
        {
            this_eq = is_zero(GR_MAT_ENTRY(mat, i, 0, sz), FLINT_MIN(i, ac), ctx);

            if (this_eq == T_FALSE)
                return T_FALSE;

            if (this_eq == T_UNKNOWN)
                eq = T_UNKNOWN;
        }

        /* main diagonal */
        if (i < ac)
        {
            this_eq = is_neg_one(GR_MAT_ENTRY(mat, i, i, sz), ctx);

            if (this_eq == T_FALSE)
                return T_FALSE;

            if (this_eq == T_UNKNOWN)
                eq = T_UNKNOWN;
        }

        /* lower triangular */
        if (i < ac - 1)
        {
            this_eq = is_zero(GR_MAT_ENTRY(mat, i, i + 1, sz), ac - i - 1, ctx);

            if (this_eq == T_FALSE)
                return T_FALSE;

            if (this_eq == T_UNKNOWN)
                eq = T_UNKNOWN;
        }
    }

    return eq;
}
