/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"


static len_t
_bell_series_cutoff(len_t n)
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
_mpz_bell_bsplit(mpz_t P, mpz_t Q, len_t a, len_t b, len_t n, len_t bmax)
{
    if (b - a < 20)
    {
        mpz_t u;
        len_t k;
        mpz_init(u);
        mpz_set_ui(P, 0UL);
        mpz_set_ui(Q, 0UL);
        mpz_set_ui(Q, (b - 1 == bmax) ? 1UL : b);
        for (k = b - 1; k >= a; k--)
        {
            mpz_set_ui(u, k);
            mpz_pow_ui(u, u, n);
            mpz_addmul(P, Q, u);
            if (k != a)
                mpz_mul_ui(Q, Q, k);
        }
        mpz_clear(u);
    }
    else
    {
        len_t m;
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
    len_t N, prec;
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
