/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_integral(arb_ptr res, arb_srcptr poly, slong len, slong prec)
{
    slong k = len - 1;

    for (k = len - 1; k > 0; k--)
        arb_div_ui(res + k, poly + k - 1, k, prec);

    arb_zero(res);
}

void
arb_poly_integral(arb_poly_t res, const arb_poly_t poly, slong prec)
{
    arb_poly_fit_length(res, poly->length + 1);
    _arb_poly_integral(res->coeffs, poly->coeffs, poly->length + 1, prec);
    _arb_poly_set_length(res, poly->length + 1);
    _arb_poly_normalise(res);
}
