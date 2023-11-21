/*
    Copyright (C) 2023 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "gr_poly.h"

void
_acb_poly_nth_derivative(acb_ptr res, acb_srcptr poly, ulong n, slong len, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, prec);
    if (_gr_poly_nth_derivative(res, poly, n, len, ctx) != GR_SUCCESS)
        _acb_vec_indeterminate(res, n);
}

void
acb_poly_nth_derivative(acb_poly_t res, const acb_poly_t poly, ulong n, slong prec)
{
    const slong len = poly->length;

    if (len <= n)
    {
        acb_poly_zero(res);
    }
    else if (n == 0)
    {
        acb_poly_set(res, poly);
    }
    else if (n == 1)
    {
        acb_poly_fit_length(res, len - 1);
        _acb_poly_derivative(res->coeffs, poly->coeffs, len, prec);
        _acb_poly_set_length(res, len - 1);
    }
    else
    {
        acb_poly_fit_length(res, len - n);
        _acb_poly_nth_derivative(res->coeffs, poly->coeffs, n, len, prec);
        _acb_poly_set_length(res, len - n);
    }
}
