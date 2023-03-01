/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
acb_poly_set_coeff_acb(acb_poly_t poly, slong n, const acb_t x)
{
    acb_poly_fit_length(poly, n + 1);

    if (n + 1 > poly->length)
    {
        _acb_vec_zero(poly->coeffs + poly->length, n - poly->length);
        poly->length = n + 1;
    }

    acb_set(poly->coeffs + n, x);
    _acb_poly_normalise(poly);
}

