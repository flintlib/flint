/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "gr_vec.h"
#include "gr_poly.h"

int
gr_poly_mul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0 || gr_is_zero(c, ctx) == T_TRUE)
        return gr_poly_zero(res, ctx);

    if (res != poly)
    {
        gr_poly_fit_length(res, len, ctx);
        _gr_poly_set_length(res, len, ctx);
    }

    status = _gr_vec_mul_scalar(res->coeffs, poly->coeffs, len, c, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_scalar_mul(gr_poly_t res, gr_srcptr c, const gr_poly_t poly, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0 || gr_is_zero(c, ctx) == T_TRUE)
        return gr_poly_zero(res, ctx);

    if (res != poly)
    {
        gr_poly_fit_length(res, len, ctx);
        _gr_poly_set_length(res, len, ctx);
    }

    status = _gr_scalar_mul_vec(res->coeffs, c, poly->coeffs, len, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_mul_ui(gr_poly_t res, const gr_poly_t poly, ulong c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0 || c == 0)
        return gr_poly_zero(res, ctx);

    if (res != poly)
    {
        gr_poly_fit_length(res, len, ctx);
        _gr_poly_set_length(res, len, ctx);
    }

    status = _gr_vec_mul_scalar_ui(res->coeffs, poly->coeffs, len, c, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_mul_si(gr_poly_t res, const gr_poly_t poly, slong c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0 || c == 0)
        return gr_poly_zero(res, ctx);

    if (res != poly)
    {
        gr_poly_fit_length(res, len, ctx);
        _gr_poly_set_length(res, len, ctx);
    }

    status = _gr_vec_mul_scalar_si(res->coeffs, poly->coeffs, len, c, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_mul_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0 || fmpz_is_zero(c))
        return gr_poly_zero(res, ctx);

    if (res != poly)
    {
        gr_poly_fit_length(res, len, ctx);
        _gr_poly_set_length(res, len, ctx);
    }

    status = _gr_vec_mul_scalar_fmpz(res->coeffs, poly->coeffs, len, c, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_mul_fmpq(gr_poly_t res, const gr_poly_t poly, const fmpq_t c, gr_ctx_t ctx)
{
    int status;
    slong len = poly->length;

    if (len == 0 || fmpq_is_zero(c))
        return gr_poly_zero(res, ctx);

    if (res != poly)
    {
        gr_poly_fit_length(res, len, ctx);
        _gr_poly_set_length(res, len, ctx);
    }

    status = _gr_vec_mul_scalar_fmpq(res->coeffs, poly->coeffs, len, c, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
