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

/* Defined in t-factor.c and t-factor_trial.c */
#define check check_factor
void check(fmpz_t n)
{
    fmpz_factor_t factor;
    fmpz_t m;
    slong i, j;

    fmpz_factor_init(factor);
    fmpz_init(m);

    fmpz_factor(factor, n);
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

    for (i = 0; i < factor->num; i++)
        for (j = i + 1; j < factor->num; j++)
        {
            if (!fmpz_cmp(factor->p + i, factor->p + j))
            {
                flint_printf("ERROR: duplicated prime factors, the form is not canonical!\n");

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

TEST_FUNCTION_START(fmpz_factor, state)
{
    int i, j, k;
    fmpz_t x, y, z, n;
    fmpz_factor_t factors;
    mpz_t y1;

    fmpz_init(x);
    mpz_init(y1);

    /* Fredrik's example */
#if FLINT64
    fmpz_set_ui(x, 3);
    fmpz_mul_ui(x, x, 73);
    fmpz_mul_ui(x, x, 137);
    fmpz_mul_ui(x, x, 1676321);
    fmpz_mul_ui(x, x, 1676321);
    fmpz_mul_ui(x, x, 1676321);
    fmpz_mul_ui(x, x, 5964848081);
    fmpz_mul_ui(x, x, 78875943472201);
    fmpz_mul_ui(x, x, 78875943472201);
    check(x);
#endif

    /* Some corner cases */
    fmpz_set_ui(x, UWORD_MAX);
    check(x);
    fmpz_set_si(x, WORD_MAX);
    check(x);
    fmpz_set_si(x, WORD_MIN);
    check(x);
    fmpz_set_si(x, COEFF_MAX);
    check(x);
    fmpz_set_si(x, COEFF_MIN);
    check(x);

    /* Small integers */
    for (i = -10000; i < 10000; i++)
    {
        fmpz_set_si(x, i);
        check(x);
    }

    /* Powers */
    for (i = 1; i < 250; i++)
    {
        for (j = 0; j < 250; j++)
        {
            fmpz_set_ui(x, i);
            fmpz_pow_ui(x, x, j);
            check(x);
        }
    }

    /* Factorials */
    for (i = 0; i < 1000; i++)
    {
        flint_mpz_fac_ui(y1, i);
        fmpz_set_mpz(x, y1);
        check(x);
    }

    /* Powers of factorials */
    for (i = 0; i < 100; i++)
    {
        for (j = 1; j < 5; j++)
        {
            flint_mpz_fac_ui(y1, i);
            fmpz_set_mpz(x, y1);
            fmpz_pow_ui(x, x, j);
            check(x);
        }
    }

    /* Whole limbs */
    for (i = 0; i < 1000; i++)
    {
        fmpz_set_ui(x, n_randtest(state));
        if (n_randint(state, 2))
            fmpz_neg(x, x);
        check(x);
    }

    /* Large negative integers */
    fmpz_set_ui(x, 10);
    fmpz_pow_ui(x, x, 100);
    fmpz_neg(x, x);
    check(x);
    flint_mpz_fac_ui(y1, 50);
    mpz_neg(y1, y1);
    fmpz_set_mpz(x, y1);
    check(x);

    mpz_clear(y1);

    fmpz_init(y);
    fmpz_init(z);
    fmpz_init(n);

    for (i = 0; i < 20; i++) /* Test random n, two factors */
    {
       randprime(x, state, 50);
       randprime(y, state, 50);

       fmpz_mul(n, x, y);

       fmpz_factor_init(factors);

       fmpz_factor(factors, n);

       if (factors->num < 2)
       {
          flint_printf("FAIL:\n");
          flint_printf("%ld factors found\n", factors->num);
          fflush(stdout);
          flint_abort();
       }

       fmpz_factor_clear(factors);
    }

    for (i = 0; i < 20; i++) /* Test random n, three factors */
    {
       randprime(x, state, 40);
       randprime(y, state, 40);
       randprime(z, state, 40);

       fmpz_mul(n, x, y);
       fmpz_mul(n, n, z);

       fmpz_factor_init(factors);

       fmpz_factor(factors, n);

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

       fmpz_factor(factors, n);

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

       fmpz_factor(factors, n);

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

       fmpz_factor(factors, n);

       if (factors->num < 1)
       {
          flint_printf("FAIL:\n");
          flint_printf("%ld factors found\n", factors->num);
          fflush(stdout);
          flint_abort();
       }
       fmpz_factor_clear(factors);
    }

    for (i = 0; i < 5; i++) /* Test random p1*p2*p3^2 */
    {
       randprime(x, state, 40);
       randprime(y, state, 40);
       randprime(z, state, 40);

       fmpz_mul(n, x, y);
       fmpz_mul(n, n, z);
       fmpz_mul(n, n, z);

       fmpz_factor_init(factors);

       fmpz_factor(factors, n);

       if (factors->num != 3)
       {
          flint_printf("FAIL:\n");
          flint_printf("%ld factors found\n", factors->num);
          fflush(stdout);
          flint_abort();
       }
       fmpz_factor_clear(factors);
    }

    for (i = 0; i < 15; i++) /* p1^e1 * p2^e2 * p3^e3 * p4^e4, e1, .., e4 in [1, .., 5] */
    {
        fmpz_set_ui(n, 1);
        for (j = 0; j < 4; j++)
        {
            slong exp;
            randprime(x, state, 20 + n_randint(state, 30));
            exp = n_randint(state, 5) + 1;
            for (k = 0; k < exp; k++)
                fmpz_mul(n, n, x);
        }

        fmpz_factor_init(factors);

        fmpz_factor(factors, n);

        if (factors->num != 4)
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
#undef check
