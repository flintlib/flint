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
gr_poly_set_coeff_scalar(gr_poly_t poly, slong n, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    gr_poly_fit_length(poly, n + 1, ctx);

    if (n + 1 > poly->length)
    {
        status |= _gr_vec_zero(GR_ENTRY(poly->coeffs, poly->length, sz), n - poly->length, ctx);
        poly->length = n + 1;
    }

    status |= gr_set(GR_ENTRY(poly->coeffs, n, sz), x, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}

int
gr_poly_set_coeff_si(gr_poly_t poly, slong n, slong x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    gr_poly_fit_length(poly, n + 1, ctx);

    if (n + 1 > poly->length)
    {
        status |= _gr_vec_zero(GR_ENTRY(poly->coeffs, poly->length, sz), n - poly->length, ctx);
        poly->length = n + 1;
    }

    status |= gr_set_si(GR_ENTRY(poly->coeffs, n, sz), x, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}

int
gr_poly_set_coeff_ui(gr_poly_t poly, slong n, ulong x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    gr_poly_fit_length(poly, n + 1, ctx);

    if (n + 1 > poly->length)
    {
        status |= _gr_vec_zero(GR_ENTRY(poly->coeffs, poly->length, sz), n - poly->length, ctx);
        poly->length = n + 1;
    }

    status |= gr_set_ui(GR_ENTRY(poly->coeffs, n, sz), x, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}

int
gr_poly_set_coeff_fmpz(gr_poly_t poly, slong n, const fmpz_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    gr_poly_fit_length(poly, n + 1, ctx);

    if (n + 1 > poly->length)
    {
        status |= _gr_vec_zero(GR_ENTRY(poly->coeffs, poly->length, sz), n - poly->length, ctx);
        poly->length = n + 1;
    }

    status |= gr_set_fmpz(GR_ENTRY(poly->coeffs, n, sz), x, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}

int
gr_poly_set_coeff_fmpq(gr_poly_t poly, slong n, const fmpq_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    gr_poly_fit_length(poly, n + 1, ctx);

    if (n + 1 > poly->length)
    {
        status |= _gr_vec_zero(GR_ENTRY(poly->coeffs, poly->length, sz), n - poly->length, ctx);
        poly->length = n + 1;
    }

    status |= gr_set_fmpq(GR_ENTRY(poly->coeffs, n, sz), x, ctx);
    _gr_poly_normalise(poly, ctx);
    return status;
}
