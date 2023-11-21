/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_derivative(acb_ptr res, acb_srcptr poly, slong len, slong prec)
{
    slong i;

    for (i = 1; i < len; i++)
        acb_mul_ui(res + i - 1, poly + i, i, prec);
}

void
acb_poly_derivative(acb_poly_t res, const acb_poly_t poly, slong prec)
{
    slong len = poly->length;

    if (len < 2)
    {
        acb_poly_zero(res);
    }
    else
    {
        acb_poly_fit_length(res, len - 1);
        _acb_poly_derivative(res->coeffs, poly->coeffs, len, prec);
        _acb_poly_set_length(res, len - 1);
    }
}
