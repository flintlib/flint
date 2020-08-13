/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void arith_euler_polynomial(fmpq_poly_t poly, ulong n)
{
    fmpz_t t;
    slong k;

    if (n == 0)
    {
        fmpq_poly_set_ui(poly, UWORD(1));
        return;
    }

    arith_bernoulli_polynomial(poly, n + 1);

    fmpz_init(t);
    fmpz_set_si(t, WORD(-2));
    for (k = n; k >= 0; k--)
    {
        fmpz_mul(poly->coeffs + k, poly->coeffs + k, t);
        fmpz_mul_ui(t, t, UWORD(2));
        fmpz_sub_ui(t, t, UWORD(2));
    }
    fmpz_zero(poly->coeffs + n + 1);
    fmpz_mul_ui(poly->den, poly->den, n + 1);
    fmpq_poly_canonicalise(poly);
    fmpz_clear(t);
}
