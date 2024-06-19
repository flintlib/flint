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

TEST_FUNCTION_START(fmpz_factor_smooth, state)
{
    int i;
    fmpz_t x, y, z, n;
    fmpz_factor_t factors;

    fmpz_init(x);

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
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong b = n_randint(state, 10000);
        if (n_randint(state, 2))
            b = -b;

        fmpz_set_si(x, b);
        checkb(x, 16);
    }

    /* Powers */
    for (i = 1; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_set_ui(x, n_randint(state, 250));
        fmpz_pow_ui(x, x, n_randint(state, 250));
        checkb(x, 10);
    }

    /* Factorials */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_fac_ui(x, n_randint(state, 1000));
        checkb(x, 12);
    }

    /* Powers of factorials */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_fac_ui(x, n_randint(state, 100));
        fmpz_pow_ui(x, x, n_randint(state, 5));
        checkb(x, 12);
    }

    /* Whole limbs */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_set_ui(x, n_randtest(state));
        if (n_randint(state, 2))
            fmpz_neg(x, x);
        checkb(x, 32);
    }

    /* Large negative integers */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_set_ui(x, 10 + n_randint(state, 10));
        fmpz_pow_ui(x, x, n_randint(state, 100));
        fmpz_neg(x, x);
        checkb(x, 8);
    }

    fmpz_init(y);
    fmpz_init(z);
    fmpz_init(n);

    for (i = 0; i < 20; i++) /* Test random n, three factors */
    {
       fmpz_randprime(x, state, 40, 0);
       fmpz_randprime(y, state, 40, 0);
       fmpz_randprime(z, state, 40, 0);

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
       fmpz_randprime(x, state, 10, 0);
       fmpz_randprime(y, state, 10, 0);
       fmpz_randprime(z, state, 40, 0);

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
       fmpz_randprime(x, state, 40, 0);

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
       fmpz_randprime(x, state, 40, 0);

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
