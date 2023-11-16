/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "qqbar.h"

void
qqbar_neg(qqbar_t res, const qqbar_t x)
{
    slong i;
    fmpz_poly_set(QQBAR_POLY(res), QQBAR_POLY(x));

    for (i = fmpz_poly_degree(QQBAR_POLY(res)) - 1; i >= 0; i -= 2)
        fmpz_neg(QQBAR_COEFFS(res) + i, QQBAR_COEFFS(res) + i);

    acb_neg(QQBAR_ENCLOSURE(res), QQBAR_ENCLOSURE(x));
}

