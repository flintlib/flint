/*
    Copyright (C) 2023 David Berghaus

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

#include <stdio.h>

int
_gr_poly_evaluate_modular(gr_ptr y, gr_srcptr poly,
    slong len, gr_srcptr x, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

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
        slong j, k, coeff_index, l, m;
        gr_ptr xs, ys, tmp, partial_results;
        gr_ptr tmp_gr;
        k = n_sqrt(len)+1;
        j = (len + k - 1) / k;

        GR_TMP_INIT(tmp_gr, ctx);
        GR_TMP_INIT_VEC(xs, j, ctx);
        GR_TMP_INIT_VEC(ys, k, ctx);
        GR_TMP_INIT_VEC(tmp, k, ctx);
        GR_TMP_INIT_VEC(partial_results, j, ctx);

        status |= _gr_vec_set_powers(xs, x, j, ctx);
        status |= gr_mul(tmp_gr, x, GR_ENTRY(xs, j-1, sz), ctx); //This is x^j
        status |= _gr_vec_set_powers(ys, tmp_gr, k, ctx);

        for (l = 0; l < j; l++){
            for (m = 0; m < k; m++){
                coeff_index = j*m+l;
                if (coeff_index < len){ //Otherwise coeffs are left to be zero
                    gr_swap(GR_ENTRY(tmp, m, sz), GR_ENTRY(poly, coeff_index, sz), ctx);
                }
            }
            status |= _gr_vec_dot(GR_ENTRY(partial_results, l, sz), NULL, 0, tmp, ys, k, ctx);
            for (m = 0; m < k; m++){ //Swap back
                coeff_index = j*m+l;
                if (coeff_index < len){
                    gr_swap(GR_ENTRY(tmp, m, sz), GR_ENTRY(poly, coeff_index, sz), ctx);
                }
            }
        }
        status |=_gr_vec_dot(y, NULL, 0, partial_results, xs, j, ctx);
        GR_TMP_CLEAR(tmp_gr, ctx);
        GR_TMP_CLEAR_VEC(xs, j, ctx);
        GR_TMP_CLEAR_VEC(ys, k, ctx);
        GR_TMP_CLEAR_VEC(tmp, k, ctx);
        GR_TMP_CLEAR_VEC(partial_results, j, ctx);
    }
    return status;
}


int
gr_poly_evaluate_modular(gr_ptr res, const gr_poly_t f, gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_modular(res, f->coeffs, f->length, x, ctx);
}
