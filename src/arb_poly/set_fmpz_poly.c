/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "arb_poly.h"

void
arb_poly_set_fmpz_poly(arb_poly_t poly, const fmpz_poly_t src, slong prec)
{
    slong i, len = fmpz_poly_length(src);

    arb_poly_fit_length(poly, len);
    _arb_poly_set_length(poly, len);

    for (i = 0; i < len; i++)
        arb_set_round_fmpz(poly->coeffs + i, src->coeffs + i, prec);
}
