/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_tan_series_basecase(gr_ptr g, gr_srcptr h, slong hlen, slong len, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    hlen = FLINT_MIN(hlen, len);

    if (hlen == 1)
    {
        status |= gr_tan(g, h, ctx);
        status |= _gr_vec_zero(GR_ENTRY(g, 1, sz), len - 1, ctx);
    }
    else if (len == 2)
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status |= gr_tan(g, h, ctx);
        status |= gr_mul(t, g, g, ctx);
        status |= gr_add_ui(t, t, 1, ctx);
        status |= gr_mul(GR_ENTRY(g, 1, sz), t, GR_ENTRY(h, 1, sz), ctx);  /* safe since hlen >= 2 */
        GR_TMP_CLEAR(t, ctx);
    }
    else
    {
        gr_ptr t, u;

        GR_TMP_INIT_VEC(t, 2 * len, ctx);
        u = GR_ENTRY(t, len, sz);

        status |= _gr_poly_sin_cos_series_basecase(t, u, h, hlen, len, 0, ctx);
        status |= _gr_poly_div_series(g, t, len, u, len, len, ctx);

        GR_TMP_CLEAR_VEC(t, 2 * len, ctx);
    }

    return status;
}

int
gr_poly_tan_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_tan_series_basecase(res->coeffs, f->coeffs, flen, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
