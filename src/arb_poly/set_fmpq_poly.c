/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_poly.h"
#include "arb_poly.h"

void
arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src, slong prec)
{
    slong i, len = src->length;

    arb_poly_fit_length(poly, len);
    _arb_poly_set_length(poly, len);

    for (i = 0; i < len; i++)
    {
        fmpq t;

        t.num = (fmpz) src->coeffs[i];
        t.den = (fmpz) *src->den;

        arb_set_fmpq(poly->coeffs + i, &t, prec);
    }
}
