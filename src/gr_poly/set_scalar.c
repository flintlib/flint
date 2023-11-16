/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "gr_poly.h"

int
gr_poly_set_scalar(gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx)
{
    if (gr_is_zero(x, ctx) == T_TRUE)
    {
        return gr_poly_zero(poly, ctx);
    }
    else
    {
        int status;
        gr_poly_fit_length(poly, 1, ctx);
        status = gr_set(poly->coeffs, x, ctx);
        _gr_poly_set_length(poly, 1, ctx);
        return status;
    }
}

int
gr_poly_set_si(gr_poly_t poly, slong x, gr_ctx_t ctx)
{
    if (x == 0)
    {
        return gr_poly_zero(poly, ctx);
    }
    else
    {
        int status;
        gr_poly_fit_length(poly, 1, ctx);
        status = gr_set_si(poly->coeffs, x, ctx);
        _gr_poly_set_length(poly, 1, ctx);
        _gr_poly_normalise(poly, ctx);
        return status;
    }
}

int
gr_poly_set_ui(gr_poly_t poly, ulong x, gr_ctx_t ctx)
{
    if (x == 0)
    {
        return gr_poly_zero(poly, ctx);
    }
    else
    {
        int status;
        gr_poly_fit_length(poly, 1, ctx);
        status = gr_set_ui(poly->coeffs, x, ctx);
        _gr_poly_set_length(poly, 1, ctx);
        _gr_poly_normalise(poly, ctx);
        return status;
    }
}

int
gr_poly_set_fmpz(gr_poly_t poly, const fmpz_t x, gr_ctx_t ctx)
{
    if (fmpz_is_zero(x))
    {
        return gr_poly_zero(poly, ctx);
    }
    else
    {
        int status;
        gr_poly_fit_length(poly, 1, ctx);
        status = gr_set_fmpz(poly->coeffs, x, ctx);
        _gr_poly_set_length(poly, 1, ctx);
        _gr_poly_normalise(poly, ctx);
        return status;
    }
}

int
gr_poly_set_fmpq(gr_poly_t poly, const fmpq_t x, gr_ctx_t ctx)
{
    if (fmpq_is_zero(x))
    {
        return gr_poly_zero(poly, ctx);
    }
    else
    {
        int status;
        gr_poly_fit_length(poly, 1, ctx);
        status = gr_set_fmpq(poly->coeffs, x, ctx);
        _gr_poly_set_length(poly, 1, ctx);
        _gr_poly_normalise(poly, ctx);
        return status;
    }
}
