/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
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

TEST_FUNCTION_START(fmpz_factor, state)
{
    int i;
    fmpz_t x, y, z, n;
    fmpz_factor_t factors;

    fmpz_init(x);

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
    fmpz_set_si(x, 1);
    check(x);
    fmpz_set_si(x, 0);
    check(x);
    fmpz_set_si(x, -1);
    check(x);

    /* Small integers */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong b = n_randint(state, 10000);
        if (n_randint(state, 2))
            b = -b;

        fmpz_set_si(x, b);
        check(x);
    }

    /* Powers */
    for (i = 1; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_set_ui(x, n_randint(state, 300));
        fmpz_pow_ui(x, x, n_randint(state, 300));
        check(x);
    }

    /* Factorials */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_fac_ui(x, n_randint(state, 1000));
        check(x);
    }

    /* Powers of factorials */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_fac_ui(x, n_randint(state, 100));
        fmpz_pow_ui(x, x, n_randint(state, 5));
        check(x);
    }

    /* Whole limbs */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_set_ui(x, n_randtest(state));
        if (n_randint(state, 2))
            fmpz_neg(x, x);
        check(x);
    }

    /* Large negative integers */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_set_ui(x, 10 + n_randint(state, 10));
        fmpz_pow_ui(x, x, n_randint(state, 100));
        fmpz_neg(x, x);
        check(x);
    }

    fmpz_init(y);
    fmpz_init(z);
    fmpz_init(n);

    for (i = 0; i < 10 + flint_test_multiplier(); i++) /* Test random n, two factors */
    {
       fmpz_randprime(x, state, 50, 0);
       fmpz_randprime(y, state, 50, 0);

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

    for (i = 0; i < 10 + flint_test_multiplier(); i++) /* Test random n, three factors */
    {
       fmpz_randprime(x, state, 40, 0);
       fmpz_randprime(y, state, 40, 0);
       fmpz_randprime(z, state, 40, 0);

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

    for (i = 0; i < 10 + flint_test_multiplier(); i++) /* Test random n, small factors */
    {
       fmpz_randprime(x, state, 10, 0);
       fmpz_randprime(y, state, 10, 0);
       fmpz_randprime(z, state, 40, 0);

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

    for (i = 0; i < 5 + flint_test_multiplier(); i++) /* Test random squares */
    {
       fmpz_randprime(x, state, 40, 0);

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

    for (i = 0; i < 5 + flint_test_multiplier(); i++) /* Test random cubes */
    {
       fmpz_randprime(x, state, 40, 0);

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

    for (i = 0; i < 5 + flint_test_multiplier(); i++) /* Test random p1*p2*p3^2 */
    {
       fmpz_randprime(x, state, 40, 0);
       fmpz_randprime(y, state, 40, 0);
       fmpz_randprime(z, state, 40, 0);

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

    for (i = 0; i < 10 + flint_test_multiplier(); i++) /* p1^e1 * p2^e2 * p3^e3 * p4^e4, e1, .., e4 in [1, .., 5] */
    {
        slong j;

        fmpz_set_ui(n, 1);

        for (j = 0; j < 4; j++)
        {
            fmpz_randprime(x, state, 20 + n_randint(state, 30), 0);
            fmpz_pow_ui(x, x, n_randint(state, 5) + 1);
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
