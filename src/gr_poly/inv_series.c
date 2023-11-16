/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_inv_series_generic(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, gr_ctx_t ctx)
{
    /* todo */
    if (Qlen <= 8 || ctx->methods[GR_METHOD_POLY_MULLOW] == (gr_funcptr) _gr_poly_mullow_generic)
    {
        return _gr_poly_inv_series_basecase(Qinv, Q, Qlen, len, ctx);
    }
    else
    {
        return _gr_poly_inv_series_newton(Qinv, Q, Qlen, len, FLINT_MIN(10, len / 2), ctx);
    }
}

int
gr_poly_inv_series(gr_poly_t Qinv, const gr_poly_t Q, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong Qlen;

    if (len == 0)
        return gr_poly_zero(Qinv, ctx);

    Qlen = Q->length;

    if (Qlen == 0)
        return GR_DOMAIN;

    if (Qlen == 1)
        len = 1;

    if (Qinv == Q)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_inv_series(t, Q, len, ctx);
        gr_poly_swap(Qinv, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Qinv, len, ctx);
    status |= _gr_poly_inv_series(Qinv->coeffs, Q->coeffs, Q->length, len, ctx);
    _gr_poly_set_length(Qinv, len, ctx);
    _gr_poly_normalise(Qinv, ctx);
    return status;
}
