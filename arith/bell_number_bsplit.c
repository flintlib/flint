/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "arith.h"


static slong
_bell_series_cutoff(slong n)
{
    double N, log_N, log_pow, log_fac;

    N = n;
    log_N = (N==0 ? 0 : log(N));
    log_pow = n * log_N;
    log_fac = N*log_N - N;
    while (log_pow - log_fac >= -2)
    {
        N++;
        log_N = log(N);
        log_pow = n * log_N;
        log_fac += log_N;
    }
    return N;
}

static void
_mpz_bell_bsplit(mpz_t P, mpz_t Q, slong a, slong b, slong n, slong bmax)
{
    if (b - a < 20)
    {
        mpz_t u;
        slong k;
        mpz_init(u);
        flint_mpz_set_ui(P, UWORD(0));
        flint_mpz_set_ui(Q, UWORD(0));
        flint_mpz_set_ui(Q, (b - 1 == bmax) ? UWORD(1) : b);
        for (k = b - 1; k >= a; k--)
        {
            flint_mpz_set_ui(u, k);
            flint_mpz_pow_ui(u, u, n);
            mpz_addmul(P, Q, u);
            if (k != a)
                flint_mpz_mul_ui(Q, Q, k);
        }
        mpz_clear(u);
    }
    else
    {
        slong m;
        mpz_t P1, Q2;
        m = (a + b) / 2;
        mpz_init(P1);
        mpz_init(Q2);
        _mpz_bell_bsplit(P1, Q, a, m, n, bmax);
        _mpz_bell_bsplit(P, Q2, m, b, n, bmax);
        mpz_mul(Q, Q, Q2);
        mpz_addmul(P, P1, Q2);
        mpz_clear(P1);
        mpz_clear(Q2);
    }
}

void
arith_bell_number_bsplit(fmpz_t b, ulong n)
{
    slong N;
    flint_bitcnt_t prec;
    mpz_t P, Q;
    mpfr_t Pf, Qf, E, one;

    N = _bell_series_cutoff(n);

    mpz_init(P);
    mpz_init(Q);

    _mpz_bell_bsplit(P, Q, 1, N + 1, n, N);

    prec = mpz_sizeinbase(P, 2) - mpz_sizeinbase(Q, 2) + 10;

    mpfr_init2(Pf, prec);
    mpfr_init2(Qf, prec);
    mpfr_init2(E, prec);
    mpfr_init2(one, 2);

    mpfr_set_z(Pf, P, GMP_RNDN);
    mpfr_set_z(Qf, Q, GMP_RNDN);
    mpfr_set_ui(one, 1, GMP_RNDN);
    mpfr_exp(E, one, GMP_RNDN);
    mpfr_mul(Qf, Qf, E, GMP_RNDN);
    mpfr_div(Pf, Pf, Qf, GMP_RNDN);
    mpfr_get_z(P, Pf, GMP_RNDN);

    fmpz_set_mpz(b, P);

    mpfr_clear(one);
    mpfr_clear(Pf);
    mpfr_clear(Qf);
    mpfr_clear(E);
    mpz_clear(P);
    mpz_clear(Q);
}
