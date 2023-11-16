/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_special.h"
#include "gr_poly.h"

int
_gr_poly_log1p_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    gr_ptr a;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len == 0)
        return GR_SUCCESS;

    flen = FLINT_MIN(flen, len);

    GR_TMP_INIT(a, ctx);
    status |= gr_log1p(a, f, ctx);

    if (flen == 1)
    {
        status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), len - 1, ctx);
    }
    else if (len == 2)
    {
        status |= gr_add_ui(res, f, 1, ctx);
        status |= gr_div(GR_ENTRY(res, 1, sz), GR_ENTRY(f, 1, sz), res, ctx);
    }
    else if (flen == 2 || _gr_vec_is_zero(GR_ENTRY(f, 1, sz), flen - 2, ctx) == T_TRUE)  /* f = a + bx^d */
    {
        slong i, j, d = flen - 1;

        status |= gr_add_ui(res, f, 1, ctx);

        for (i = 1, j = d; j < len; j += d, i++)
        {
            if (i == 1)
                status |= gr_div(GR_ENTRY(res, j, sz), GR_ENTRY(f, d, sz), res, ctx);
            else
                status |= gr_mul(GR_ENTRY(res, j, sz), GR_ENTRY(res, j - d, sz), GR_ENTRY(res, d, sz), ctx);

            status |= _gr_vec_zero(GR_ENTRY(res, j - d + 1, sz), flen - 2, ctx);
        }

        status |= _gr_vec_zero(GR_ENTRY(res, j - d + 1, sz), len - (j - d + 1), ctx);

        for (i = 2, j = 2 * d; j < len; j += d, i++)
            status |= gr_div_si(GR_ENTRY(res, j, sz), GR_ENTRY(res, j, sz), i % 2 ? i : -i, ctx);
    }
    else
    {
        gr_ptr f_diff, f_inv;
        slong alloc;

        alloc = len + (flen - 1);
        GR_TMP_INIT_VEC(f_inv, alloc, ctx);

        f_diff = GR_ENTRY(f_inv, len, sz);

        if (status == GR_SUCCESS)
        {
            status |= _gr_poly_derivative(f_diff, f, flen, ctx);
            status |= gr_add_ui(res, f, 1, ctx);
            status |= _gr_vec_set(GR_ENTRY(res, 1, sz), GR_ENTRY(f, 1, sz), flen - 1, ctx);
            status |= _gr_poly_div_series(f_inv, f_diff, flen - 1, res, flen, len, ctx);
            status |= _gr_poly_integral(res, f_inv, len, ctx);
        }

        GR_TMP_CLEAR_VEC(f_inv, alloc, ctx);
    }

    gr_swap(res, a, ctx);

    GR_TMP_CLEAR(a, ctx);
    return status;
}

int
gr_poly_log1p_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_log1p_series(res->coeffs, f->coeffs, f->length, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
