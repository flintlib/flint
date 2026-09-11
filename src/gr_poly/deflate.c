/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

ulong
gr_poly_deflation(const gr_poly_t poly, gr_ctx_t ctx)
{
    slong i, sz = ctx->sizeof_elem;
    ulong deflation;

    if (poly->length <= 1)
        return poly->length;

    /* The exponents of the nonzero terms. A term whose coefficient cannot
       be verified to be zero is treated as nonzero, so the deflation is
       always valid (possibly not optimal). */
    deflation = 0;
    for (i = 1; i < poly->length - 1 && deflation != 1; i++)
        if (gr_is_zero(GR_ENTRY(poly->coeffs, i, sz), ctx) != T_TRUE)
            deflation = n_gcd(deflation, i);

    return n_gcd(deflation, poly->length - 1);
}

int
gr_poly_deflate(gr_poly_t res, const gr_poly_t poly, ulong deflation, gr_ctx_t ctx)
{
    slong i, len, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (deflation == 0)
        return GR_DOMAIN;

    if (poly->length <= 1 || deflation == 1)
        return gr_poly_set(res, poly, ctx);

    len = (poly->length - 1) / deflation + 1;

    if (res == poly)
    {
        for (i = 1; i < len; i++)
            gr_swap(GR_ENTRY(res->coeffs, i, sz), GR_ENTRY(res->coeffs, i * deflation, sz), ctx);
    }
    else
    {
        gr_poly_fit_length(res, len, ctx);
        for (i = 0; i < len; i++)
            status |= gr_set(GR_ENTRY(res->coeffs, i, sz), GR_ENTRY(poly->coeffs, i * deflation, sz), ctx);
    }

    _gr_poly_set_length_normalise(res, len, ctx);

    return status;
}

int
gr_poly_inflate(gr_poly_t res, const gr_poly_t poly, ulong inflation, gr_ctx_t ctx)
{
    slong len;
    int status = GR_SUCCESS;

    if (inflation == 0)
    {
        /* poly(x^0) = poly(1) */
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status |= gr_one(t, ctx);
        status |= gr_poly_evaluate(t, poly, t, ctx);
        status |= gr_poly_set_scalar(res, t, ctx);
        GR_TMP_CLEAR(t, ctx);
        return status;
    }

    if (poly->length <= 1 || inflation == 1)
        return gr_poly_set(res, poly, ctx);

    len = (poly->length - 1) * inflation + 1;

    if (res == poly)
    {
        gr_poly_fit_length(res, len, ctx);
        status |= _gr_poly_inflate(res->coeffs, poly->length, inflation, ctx);
    }
    else
    {
        gr_poly_fit_length(res, len, ctx);
        status |= _gr_vec_set(res->coeffs, poly->coeffs, poly->length, ctx);
        status |= _gr_poly_inflate(res->coeffs, poly->length, inflation, ctx);
    }

    _gr_poly_set_length_normalise(res, len, ctx);

    return status;
}
