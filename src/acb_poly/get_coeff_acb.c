/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
acb_poly_get_coeff_acb(acb_t x, const acb_poly_t poly, slong n)
{
    if (n < poly->length)
        acb_set(x, poly->coeffs + n);
    else
        acb_zero(x);
}

