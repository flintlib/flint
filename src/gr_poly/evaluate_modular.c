/*
    Copyright (C) 2023 David Berghaus

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_evaluate_modular(gr_ptr y, gr_srcptr poly,
    slong len, gr_srcptr x, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_method_void_unary_op set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);

    if (len <= 2)
    {
        if (len == 0)
            return gr_zero(y, ctx);

        if (len == 1)
            return gr_set(y, poly, ctx);

        status |= gr_mul(y, x, GR_ENTRY(poly, 1, sz), ctx);
        status |= gr_add(y, y, GR_ENTRY(poly, 0, sz), ctx);
    }
    else
    {
        slong i, j, k, coeff_index, l, m;
        gr_ptr xs, ys, tmp, partial_results;
        k = n_sqrt(len)+1;
        j = (len + k - 1) / k;

        TMP_INIT;

        TMP_START;
        tmp = TMP_ALLOC(sz * k);
        GR_TMP_INIT_VEC(xs, j + 1, ctx);
        GR_TMP_INIT_VEC(ys, k, ctx);
        GR_TMP_INIT_VEC(partial_results, j, ctx);

        status |= _gr_vec_set_powers(xs, x, j + 1, ctx);
        status |= _gr_vec_set_powers(ys, GR_ENTRY(xs, j, sz), k, ctx);

        for (l = 0; l < j; l++)
        {
            i = 0; /* Count number of coeffs in this row */
            for (m = 0; m < k; m++)
            {
                coeff_index = j*m+l;
                if (coeff_index < len)
                { 
                    set_shallow(GR_ENTRY(tmp, m, sz), GR_ENTRY(poly, coeff_index, sz), ctx);
                    i++;
                }
                else
                {
                    break;
                }
            }
            status |= _gr_vec_dot(GR_ENTRY(partial_results, l, sz), NULL, 0, tmp, ys, i, ctx);
        }
        status |=_gr_vec_dot(y, NULL, 0, partial_results, xs, j, ctx);
        GR_TMP_CLEAR_VEC(xs, j + 1, ctx);
        GR_TMP_CLEAR_VEC(ys, k, ctx);
        GR_TMP_CLEAR_VEC(partial_results, j, ctx);
        TMP_END;
    }
    return status;
}

int
gr_poly_evaluate_modular(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_modular(res, f->coeffs, f->length, x, ctx);
}
