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
_gr_poly_tan_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i, m, n;
    gr_ptr t, u;
    slong a[FLINT_BITS];

    flen = FLINT_MIN(flen, len);

    if (len < cutoff)
        return _gr_poly_tan_series_basecase(res, f, flen, len, ctx);

    cutoff = FLINT_MAX(cutoff, 2);

    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = (n + 1) / 2);

    status |= _gr_poly_tan_series_basecase(res, f, flen, n, ctx);

    if (status != GR_SUCCESS)
        return status;

    GR_TMP_INIT_VEC(t, 2 * len, ctx);
    u = GR_ENTRY(t, len, sz);

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];

        status |= _gr_poly_mullow(u, res, m, res, m, n, ctx);
        status |= gr_add_ui(u, u, 1, ctx);
        status |= _gr_poly_atan_series(t, res, m, n, ctx);
        status |= _gr_poly_sub(GR_ENTRY(t, m, sz), GR_ENTRY(f, m, sz), FLINT_MAX(0, flen - m), GR_ENTRY(t, m, sz), n - m, ctx);
        status |= _gr_poly_mullow(GR_ENTRY(res, m, sz), u, n, GR_ENTRY(t, m, sz), n - m, n - m, ctx);
    }

    GR_TMP_CLEAR_VEC(t, 2 * len, ctx);

    return status;
}

int
gr_poly_tan_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (flen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (flen == 1)
        len = 1;

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_tan_series_newton(res->coeffs, f->coeffs, flen, len, cutoff, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
