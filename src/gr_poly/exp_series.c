/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

/* todo */
int
_gr_poly_exp_series_generic(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
{
    hlen = FLINT_MIN(hlen, n);

    if (n < 10 || hlen < 10 || ctx->methods[GR_METHOD_POLY_MULLOW] == (gr_funcptr) _gr_poly_mullow_generic)
    {
        return _gr_poly_exp_series_basecase(f, h, hlen, n, ctx);
    }
    else if (n < 30 && hlen > 0.9 * n)
    {
        return _gr_poly_exp_series_basecase_mul(f, h, hlen, n, ctx);
    }
    else
    {
        return _gr_poly_exp_series_newton(f, NULL, h, hlen, n, 30, ctx);
    }
}

int
gr_poly_exp_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen = h->length;

    if (n == 0)
        return gr_poly_zero(f, ctx);

    if (hlen == 0)
        return gr_poly_one(f, ctx);

    gr_poly_fit_length(f, n, ctx);
    status |= _gr_poly_exp_series(f->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(f, n, ctx);
    _gr_poly_normalise(f, ctx);
    return status;
}
