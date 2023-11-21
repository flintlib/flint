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
_gr_poly_log_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len == 0)
        return GR_SUCCESS;

    flen = FLINT_MIN(flen, len);

    if (flen == 1)
    {
        status |= gr_log(res, f, ctx);
        status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), len - 1, ctx);
    }
    else if (len == 2)
    {
        status |= gr_div(GR_ENTRY(res, 1, sz), GR_ENTRY(f, 1, sz), f, ctx);
        status |= gr_log(res, f, ctx);
    }
    else if (flen == 2 || _gr_vec_is_zero(GR_ENTRY(f, 1, sz), flen - 2, ctx) == T_TRUE)  /* f = a + bx^d */
    {
        slong i, j, d = flen - 1;

        for (i = 1, j = d; j < len; j += d, i++)
        {
            if (i == 1)
                status |= gr_div(GR_ENTRY(res, j, sz), GR_ENTRY(f, d, sz), f, ctx);
            else
                status |= gr_mul(GR_ENTRY(res, j, sz), GR_ENTRY(res, j - d, sz), GR_ENTRY(res, d, sz), ctx);

            status |= _gr_vec_zero(GR_ENTRY(res, j - d + 1, sz), flen - 2, ctx);
        }

        status |= _gr_vec_zero(GR_ENTRY(res, j - d + 1, sz), len - (j - d + 1), ctx);

        for (i = 2, j = 2 * d; j < len; j += d, i++)
            status |= gr_div_si(GR_ENTRY(res, j, sz), GR_ENTRY(res, j, sz), i % 2 ? i : -i, ctx);

        status |= gr_log(res, f, ctx); /* done last to allow aliasing */
    }
    else
    {
        gr_ptr f_diff, f_inv;
        gr_ptr a;
        slong alloc;

        alloc = (len - 1) + (flen - 1) + 1;

        GR_TMP_INIT_VEC(f_inv, alloc, ctx);

        f_diff = GR_ENTRY(f_inv, len - 1, sz);

        a = GR_ENTRY(f_diff, flen - 1, sz);

        status |= gr_log(a, f, ctx);

        if (status == GR_SUCCESS)
        {
            status |= _gr_poly_derivative(f_diff, f, flen, ctx);
            status |= _gr_poly_div_series(f_inv, f_diff, flen - 1, f, flen, len - 1, ctx);
            status |= _gr_poly_integral(res, f_inv, len, ctx);
            gr_swap(res, a, ctx);
        }

        GR_TMP_CLEAR_VEC(f_inv, alloc, ctx);
    }

    return status;
}

int
gr_poly_log_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 0)
        return GR_DOMAIN;

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_log_series(res->coeffs, f->coeffs, f->length, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
