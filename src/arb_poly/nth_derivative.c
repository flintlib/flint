/*
    Copyright (C) 2023 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

void
_arb_poly_nth_derivative(arb_ptr res, arb_srcptr poly, ulong n, slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_real_arb(ctx, prec);
    if (_gr_poly_nth_derivative(res, poly, n, len, ctx) != GR_SUCCESS)
        _arb_vec_indeterminate(res, n);
}

void
arb_poly_nth_derivative(arb_poly_t res, const arb_poly_t poly, ulong n, slong prec)
{
    const slong len = poly->length;

    if (len <= n)
    {
        arb_poly_zero(res);
    }
    else if (n == 0)
    {
        arb_poly_set(res, poly);
    }
    else if (n == 1)
    {
        arb_poly_fit_length(res, len - 1);
        _arb_poly_derivative(res->coeffs, poly->coeffs, len, prec);
        _arb_poly_set_length(res, len - 1);
    }
    else
    {
        arb_poly_fit_length(res, len - n);
        _arb_poly_nth_derivative(res->coeffs, poly->coeffs, n, len, prec);
        _arb_poly_set_length(res, len - n);
    }
}
