/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_set_array_mpq(fmpz * poly, fmpz_t den, const mpq_t * a, slong n)
{
    slong i;
    mpz_t d, t;

    flint_mpz_init_set_ui(d, 1);
    mpz_init(t);
    for (i = 0; i < n; i++)
    {
        mpz_lcm(d, d, mpq_denref(a[i]));
    }

    for (i = 0; i < n; i++)
    {
        __mpz_struct *ptr = _fmpz_promote(poly + i);

        mpz_divexact(t, d, mpq_denref(a[i]));
        mpz_mul(ptr, mpq_numref(a[i]), t);
        _fmpz_demote_val(poly + i);
    }

    fmpz_set_mpz(den, d);
    mpz_clear(d);
    mpz_clear(t);
}

void fmpq_poly_set_array_mpq(fmpq_poly_t poly, const mpq_t * a, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(poly);
    }
    else
    {
        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_set_array_mpq(poly->coeffs, poly->den, a, n);
        _fmpq_poly_set_length(poly, n);
        _fmpq_poly_normalise(poly);
    }
}
