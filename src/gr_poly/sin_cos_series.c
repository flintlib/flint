/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

static int
_gr_poly_sin_cos_series_with_pi(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen,
        slong n, int times_pi, gr_ctx_t ctx)
{
    slong cutoff;

    hlen = FLINT_MIN(hlen, n);

    if (hlen < 20)
        return _gr_poly_sin_cos_series_basecase(s, c, h, hlen, n, times_pi, ctx);

    if (gr_ctx_is_finite(ctx) == T_TRUE)
    {
        /* Reasonable tuning for nmod */
        cutoff = 500;
    }
    else if (gr_ctx_has_real_prec(ctx) == T_TRUE)
    {
        /* Reasonable tuning for arb */
        GR_MUST_SUCCEED(gr_ctx_get_real_prec(&cutoff, ctx));
        cutoff = FLINT_MAX(cutoff, 1);
        cutoff = FLINT_MIN(500, 20 + 100000 / cutoff);
    }
    else
    {
        /* Reasonable tuning checked for fmpq */
        cutoff = 25;
    }

    if (hlen < cutoff)
        return _gr_poly_sin_cos_series_basecase(s, c, h, hlen, n, times_pi, ctx);
    else
        return _gr_poly_sin_cos_series_newton(s, c, h, hlen, n, cutoff, times_pi, ctx);
}

int
_gr_poly_sin_cos_series(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
{
    return _gr_poly_sin_cos_series_with_pi(s, c, h, hlen, n, 0, ctx);
}

int
_gr_poly_sin_cos_pi_series(gr_ptr s, gr_ptr c, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
{
    return _gr_poly_sin_cos_series_with_pi(s, c, h, hlen, n, 1, ctx);
}

int
_gr_poly_sin_series(gr_ptr s, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
{
    return _gr_poly_sin_cos_series_with_pi(s, NULL, h, hlen, n, 0, ctx);
}

int
_gr_poly_sin_pi_series(gr_ptr s, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
{
    return _gr_poly_sin_cos_series_with_pi(s, NULL, h, hlen, n, 1, ctx);
}

int
_gr_poly_cos_series(gr_ptr c, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
{
    return _gr_poly_sin_cos_series_with_pi(NULL, c, h, hlen, n, 0, ctx);
}

int
_gr_poly_cos_pi_series(gr_ptr c, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx)
{
    return _gr_poly_sin_cos_series_with_pi(NULL, c, h, hlen, n, 1, ctx);
}

int
gr_poly_sin_cos_series(gr_poly_t s, gr_poly_t c,
        const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    slong hlen = h->length;
    int status = GR_SUCCESS;

    if (n == 0)
    {
        status |= gr_poly_zero(s, ctx);
        status |= gr_poly_zero(c, ctx);
        return status;
    }

    if (hlen == 0)
    {
        status |= gr_poly_zero(s, ctx);
        status |= gr_poly_one(c, ctx);
        return status;
    }

    if (hlen == 1)
        n = 1;

    gr_poly_fit_length(s, n, ctx);
    gr_poly_fit_length(c, n, ctx);
    status |= _gr_poly_sin_cos_series(s->coeffs, c->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(s, n, ctx);
    _gr_poly_normalise(s, ctx);
    _gr_poly_set_length(c, n, ctx);
    _gr_poly_normalise(c, ctx);
    return status;
}

int
gr_poly_sin_cos_pi_series(gr_poly_t s, gr_poly_t c,
        const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    slong hlen = h->length;
    int status = GR_SUCCESS;

    if (n == 0)
    {
        status |= gr_poly_zero(s, ctx);
        status |= gr_poly_zero(c, ctx);
        return status;
    }

    if (hlen == 0)
    {
        status |= gr_poly_zero(s, ctx);
        status |= gr_poly_one(c, ctx);
        return status;
    }

    if (hlen == 1)
        n = 1;

    gr_poly_fit_length(s, n, ctx);
    gr_poly_fit_length(c, n, ctx);
    status |= _gr_poly_sin_cos_pi_series(s->coeffs, c->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(s, n, ctx);
    _gr_poly_normalise(s, ctx);
    _gr_poly_set_length(c, n, ctx);
    _gr_poly_normalise(c, ctx);
    return status;
}

int
gr_poly_sin_series(gr_poly_t s, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    slong hlen = h->length;
    int status = GR_SUCCESS;

    if (n == 0 || hlen == 0)
        return gr_poly_zero(s, ctx);

    if (hlen == 1)
        n = 1;

    gr_poly_fit_length(s, n, ctx);
    status |= _gr_poly_sin_series(s->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(s, n, ctx);
    _gr_poly_normalise(s, ctx);
    return status;
}

int
gr_poly_sin_pi_series(gr_poly_t s, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    slong hlen = h->length;
    int status = GR_SUCCESS;

    if (n == 0 || hlen == 0)
        return gr_poly_zero(s, ctx);

    if (hlen == 1)
        n = 1;

    gr_poly_fit_length(s, n, ctx);
    status |= _gr_poly_sin_pi_series(s->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(s, n, ctx);
    _gr_poly_normalise(s, ctx);
    return status;
}

int
gr_poly_cos_series(gr_poly_t c, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    slong hlen = h->length;
    int status = GR_SUCCESS;

    if (n == 0)
        return gr_poly_zero(c, ctx);

    if (hlen == 0)
        return gr_poly_one(c, ctx);

    gr_poly_fit_length(c, n, ctx);
    status |= _gr_poly_cos_series(c->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(c, n, ctx);
    _gr_poly_normalise(c, ctx);
    return status;
}

int
gr_poly_cos_pi_series(gr_poly_t c, const gr_poly_t h, slong n, gr_ctx_t ctx)
{
    slong hlen = h->length;
    int status = GR_SUCCESS;

    if (n == 0)
        return gr_poly_zero(c, ctx);

    if (hlen == 0)
        return gr_poly_one(c, ctx);

    gr_poly_fit_length(c, n, ctx);
    status |= _gr_poly_cos_pi_series(c->coeffs, h->coeffs, hlen, n, ctx);
    _gr_poly_set_length(c, n, ctx);
    _gr_poly_normalise(c, ctx);
    return status;
}

