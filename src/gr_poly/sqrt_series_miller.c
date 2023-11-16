/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_sqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    fmpq_t q;

    status |= gr_sqrt(res, f, ctx);
    if (status != GR_SUCCESS)
        return status;

    *fmpq_numref(q) = 1;
    *fmpq_denref(q) = 2;

    if (gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
    {
        status |= _gr_vec_reciprocals(GR_ENTRY(res, 1, sz), len - 1, ctx);

        if (status == GR_SUCCESS)
            status |= _gr_poly_pow_series_fmpq_recurrence(res, f, flen, q, len, 3, ctx);
    }
    else
    {
        status = _gr_poly_pow_series_fmpq_recurrence(res, f, flen, q, len, 1, ctx);
    }

    return status;
}

int
gr_poly_sqrt_series_miller(gr_poly_t res, const gr_poly_t h, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong hlen;

    hlen = h->length;

    if (hlen == 0 || len == 0)
        return gr_poly_zero(res, ctx);

    if (hlen == 1)
        len = 1;

    if (res == h)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_sqrt_series(t, h, len, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(res, len, ctx);
    status |= _gr_poly_sqrt_series_miller(res->coeffs, h->coeffs, h->length, len, ctx);
    _gr_poly_set_length(res, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
