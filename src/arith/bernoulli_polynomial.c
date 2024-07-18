/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "arith.h"

void arith_bernoulli_polynomial(fmpq_poly_t poly, ulong n)
{
    fmpz_t t;
    fmpz * den;
    slong k;

    if (n == 0)
    {
        fmpq_poly_set_ui(poly, UWORD(1));
        return;
    }

    fmpq_poly_fit_length(poly, n + 1);

    fmpz_init(t);
    den = _fmpz_vec_init(n + 1);

    _arith_bernoulli_number_vec(poly->coeffs, den, n + 1);

    /* Multiply the odd term by binomial(n,1) = n */
    fmpz_mul_ui(poly->coeffs + 1, poly->coeffs + 1, n);

    /* Multiply even terms by binomial coefficients */
    fmpz_one(t);
    for (k = 2; k <= n; k += 2)
    {
        fmpz_mul2_uiui(t, t, n-k+1, n-k+2);
        fmpz_divexact2_uiui(t, t, k, k-1);
        fmpz_mul(poly->coeffs + k, poly->coeffs + k, t);
    }

    /* Convert to common denominator */
    fmpz_primorial(poly->den, n + 2);
    for (k = 0; k <= n; k++)
    {
        fmpz_mul(poly->coeffs + k, poly->coeffs+k, poly->den);
        fmpz_divexact(poly->coeffs + k, poly->coeffs + k, den + k);
    }

    _fmpz_poly_reverse(poly->coeffs, poly->coeffs, n + 1, n + 1);
    _fmpq_poly_set_length(poly, n + 1);
    fmpq_poly_canonicalise(poly);

    _fmpz_vec_clear(den, n + 1);
    fmpz_clear(t);
}
