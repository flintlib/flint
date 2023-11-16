/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
arb_poly_set_coeff_arb(arb_poly_t poly, slong n, const arb_t x)
{
    arb_poly_fit_length(poly, n + 1);

    if (n + 1 > poly->length)
    {
        _arb_vec_zero(poly->coeffs + poly->length, n - poly->length);
        poly->length = n + 1;
    }

    arb_set(poly->coeffs + n, x);
    _arb_poly_normalise(poly);
}

