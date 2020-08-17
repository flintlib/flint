/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void _arith_bernoulli_number_zeta(fmpz_t num, fmpz_t den, ulong n)
{
    mpz_t r;
    mpfr_t t, u, z, pi;
    flint_bitcnt_t prec, pi_prec;

    arith_bernoulli_number_denom(den, n);

    if (n % 2)
    {
        fmpz_set_si(num, -(n == 1));
        return;
    }

    if (n < BERNOULLI_SMALL_NUMER_LIMIT)
    {
        fmpz_set_si(num, _bernoulli_numer_small[n / 2]);
        return;
    }

    prec = arith_bernoulli_number_size(n) + fmpz_bits(den);
    prec += 10 + 2*FLINT_BIT_COUNT(n);
    prec = prec * 1.001;
    pi_prec = prec;

    mpz_init(r);
    mpfr_init2(t, prec);
    mpfr_init2(u, prec);
    mpfr_init2(z, prec);
    mpfr_init2(pi, pi_prec);

    /* t = 2 * n! / (2*pi)^n */
    flint_mpz_fac_ui(r, n);
    mpfr_set_z(t, r, GMP_RNDN);
    mpfr_mul_2exp(t, t, 1, GMP_RNDN);
    mpfr_const_pi(pi, GMP_RNDN);
    mpfr_mul_2exp(pi, pi, 1, GMP_RNDN);
    mpfr_pow_ui(pi, pi, n, GMP_RNDN);
    mpfr_div(t, t, pi, GMP_RNDN);

    /* t = t / zeta(n) */
    mpfr_zeta_inv_euler_product(z, n, 0);
    mpfr_div(t, t, z, GMP_RNDN);

    /* round numerator */
    fmpz_get_mpz(r, den);
    mpfr_mul_z(t, t, r, GMP_RNDN);
    mpfr_round(t, t);
    mpfr_get_z(r, t, GMP_RNDN);
    fmpz_set_mpz(num, r);

    if (n % 4 == 0)
        fmpz_neg(num, num);

    mpz_clear(r);
    mpfr_clear(t);
    mpfr_clear(u);
    mpfr_clear(z);
    mpfr_clear(pi);
}

