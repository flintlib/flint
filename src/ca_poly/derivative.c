/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_derivative(ca_ptr res, ca_srcptr poly, slong len, ca_ctx_t ctx)
{
    slong i;

    for (i = 1; i < len; i++)
        ca_mul_ui(res + (i - 1), poly + i, i, ctx);
}

void
ca_poly_derivative(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)
{
    slong len = poly->length;

    if (len < 2)
    {
        ca_poly_zero(res, ctx);
    }
    else
    {
        ca_poly_fit_length(res, len - 1, ctx);
        _ca_poly_derivative(res->coeffs, poly->coeffs, len, ctx);
        _ca_poly_set_length(res, len - 1, ctx);
    }
}
