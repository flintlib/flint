/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz_factor.h"

void checkb(fmpz_t n, slong bits)
{
    fmpz_factor_t factor;
    fmpz_t m;
    slong i;

    fmpz_factor_init(factor);
    fmpz_init(m);

    fmpz_factor_smooth(factor, n, bits, 0);
    fmpz_factor_expand(m, factor);

    if (!fmpz_equal(n, m))
    {
        flint_printf("ERROR: factors do not unfactor to original number!\n");

        flint_printf("input: ");
        fmpz_print(n);
        flint_printf("\n");

        flint_printf("computed factors: ");
        fmpz_factor_print(factor);
        flint_printf("\n");

        flint_printf("value: ");
        fmpz_print(m);
        flint_printf("\n");

        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < factor->num; i++)
    {
        if (!fmpz_is_probabprime(factor->p + i))
        {
            flint_printf("ERROR: factor is not prime!\n");

            flint_printf("input: ");
            fmpz_print(n);
            flint_printf("\n");

            flint_printf("computed factors: ");
            fmpz_factor_print(factor);
            flint_printf("\n");

            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_clear(m);
    fmpz_factor_clear(factor);
}

/* Defined in t-factor.c and t-factor_smooth.c */
#ifndef randprime
#define randprime randprime
void randprime(fmpz_t p, flint_rand_t state, slong bits)
{
    fmpz_randbits(p, state, bits);

    if (fmpz_sgn(p) < 0)
       fmpz_neg(p, p);

    if (fmpz_is_even(p))
       fmpz_add_ui(p, p, 1);

    while (!fmpz_is_probabprime(p))
       fmpz_add_ui(p, p, 2);
}
#endif

TEST_FUNCTION_START(fmpz_factor_smooth, state)
{
    int i, j;
    fmpz_t x, y, z, n;
    fmpz_factor_t factors;
    mpz_t y1;

    fmpz_init(x);
    mpz_init(y1);

    /* Some corner cases */
    fmpz_set_ui(x, UWORD_MAX);
    checkb(x, 32);
    fmpz_set_si(x, WORD_MAX);
    checkb(x, 32);
    fmpz_set_si(x, WORD_MIN);
    checkb(x, 32);
    fmpz_set_si(x, COEFF_MAX);
    checkb(x, 32);
    fmpz_set_si(x, COEFF_MIN);
    checkb(x, 32);

    /* Small integers */
    for (i = -10000; i < 10000; i++)
    {
        fmpz_set_si(x, i);
        checkb(x, 16);
    }

    /* Powers */
    for (i = 1; i < 250; i++)
    {
        for (j = 0; j < 250; j++)
        {
            fmpz_set_ui(x, i);
            fmpz_pow_ui(x, x, j);
            checkb(x, 10);
        }
    }

    /* Factorials */
    for (i = 0; i < 1000; i++)
    {
        flint_mpz_fac_ui(y1, i);
        fmpz_set_mpz(x, y1);
        checkb(x, 12);
    }

    /* Powers of factorials */
    for (i = 0; i < 100; i++)
    {
        for (j = 1; j < 5; j++)
        {
            flint_mpz_fac_ui(y1, i);
            fmpz_set_mpz(x, y1);
            fmpz_pow_ui(x, x, j);
            checkb(x, 12);
        }
    }

    /* Whole limbs */
    for (i = 0; i < 1000; i++)
    {
        fmpz_set_ui(x, n_randtest(state));
        if (n_randint(state, 2))
            fmpz_neg(x, x);
        checkb(x, 32);
    }

    /* Large negative integers */
    fmpz_set_ui(x, 10);
    fmpz_pow_ui(x, x, 100);
    fmpz_neg(x, x);
    checkb(x, 8);
    flint_mpz_fac_ui(y1, 50);
    mpz_neg(y1, y1);
    fmpz_set_mpz(x, y1);
    checkb(x, 8);

    mpz_clear(y1);

    fmpz_init(y);
    fmpz_init(z);
    fmpz_init(n);

    for (i = 0; i < 20; i++) /* Test random n, three factors */
    {
       randprime(x, state, 40);
       randprime(y, state, 40);
       randprime(z, state, 40);

       fmpz_mul(n, x, y);
       fmpz_mul(n, n, z);

       fmpz_factor_init(factors);

       fmpz_factor_smooth(factors, n, 60, 1);

       if (factors->num < 3)
       {
          flint_printf("FAIL:\n");
          flint_printf("%ld factors found\n", factors->num);
          fflush(stdout);
          flint_abort();
       }

       fmpz_factor_clear(factors);
    }

    for (i = 0; i < 20; i++) /* Test random n, small factors */
    {
       randprime(x, state, 10);
       randprime(y, state, 10);
       randprime(z, state, 40);

       fmpz_mul(n, x, y);
       fmpz_mul(n, n, z);

       fmpz_factor_init(factors);

       fmpz_factor_smooth(factors, n, 20, 1);

       if (factors->num < 3 && !fmpz_equal(x, y))
       {
          flint_printf("FAIL:\n");
          flint_printf("%ld factors found\n", factors->num);
          fflush(stdout);
          flint_abort();
       }

       fmpz_factor_clear(factors);
    }

    for (i = 0; i < 5; i++) /* Test random squares */
    {
       randprime(x, state, 40);

       fmpz_mul(n, x, x);

       fmpz_factor_init(factors);

       fmpz_factor_smooth(factors, n, 60, 1);

       if (factors->num < 1)
       {
          flint_printf("FAIL:\n");
          flint_printf("%ld factors found\n", factors->num);
          fflush(stdout);
          flint_abort();
       }

       fmpz_factor_clear(factors);
    }

    for (i = 0; i < 5; i++) /* Test random cubes */
    {
       randprime(x, state, 40);

       fmpz_mul(n, x, x);
       fmpz_mul(n, n, x);

       fmpz_factor_init(factors);

       fmpz_factor_smooth(factors, n, 60, 1);

       if (factors->num < 1)
       {
          flint_printf("FAIL:\n");
          flint_printf("%ld factors found\n", factors->num);
          fflush(stdout);
          flint_abort();
       }
       fmpz_factor_clear(factors);
    }

    fmpz_clear(n);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(z);

    TEST_FUNCTION_END(state);
}
