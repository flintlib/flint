/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2015 William Hart
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "fmpq.h"

double fmpq_get_d(const fmpq_t a)
{
    double d;
    mpq_t z;
    flint_mpq_init_set_readonly(z, a);
    d = mpq_get_d(z);
    flint_mpq_clear_readonly(z);
    return d;
}

void fmpq_get_mpz_frac(mpz_t a, mpz_t b, fmpq_t c)
{
   fmpz_get_mpz(a, fmpq_numref(c));
   fmpz_get_mpz(b, fmpq_denref(c));
}

void fmpq_get_mpq(mpq_t dest, const fmpq_t src)
{
    fmpz_get_mpz(mpq_numref(dest), fmpq_numref(src));
    fmpz_get_mpz(mpq_denref(dest), fmpq_denref(src));
}

int
fmpq_get_mpfr(mpfr_t r, const fmpq_t x, mpfr_rnd_t rnd)
{
    __mpq_struct mpq;
    fmpz p, q;
    ulong pp, qq;

    p = *fmpq_numref(x);
    q = *fmpq_denref(x);

    if (p == 0)
        return mpfr_set_ui(r, 0, rnd);

    if (COEFF_IS_MPZ(p))
        mpq._mp_num = *COEFF_TO_PTR(p);
    else
    {
        pp = FLINT_ABS(p);
        mpq._mp_num._mp_alloc = 1;
        mpq._mp_num._mp_size = (p < 0) ? -1 : 1;
        mpq._mp_num._mp_d = &pp;
    }

    if (COEFF_IS_MPZ(q))
        mpq._mp_den = *COEFF_TO_PTR(q);
    else
    {
        qq = q;
        mpq._mp_den._mp_alloc = 1;
        mpq._mp_den._mp_size = 1;
        mpq._mp_den._mp_d = &qq;
    }

    return mpfr_set_q(r, &mpq, rnd);
}
