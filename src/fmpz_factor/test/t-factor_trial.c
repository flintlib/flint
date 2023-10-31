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
#define check check_factor_trial
void check(fmpz_t n)
{
    fmpz_factor_t factor, factor2;
    fmpz_t m;
    slong i;
    int full;

    fmpz_factor_init(factor);
    fmpz_factor_init(factor2);
    fmpz_init(m);

    full = fmpz_factor_trial(factor, n, 1000);
    fmpz_factor_trial_range(factor2, n, 0, 1000);
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

    if (factor->num != factor2->num)
    {
        flint_printf("ERROR: number of factors do not agree\n");

	    flint_printf("n = ");
	    fmpz_print(n);
	    flint_printf("\n");

        flint_printf("factor_trial computed factors: ");
	    fmpz_factor_print(factor);
	    flint_printf("\n");

        flint_printf("factor_trial_range computed factors: ");
        fmpz_factor_print(factor2);
        flint_printf("\n");

        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < factor->num - (full != 1); i++)
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
    fmpz_factor_clear(factor2);
}

TEST_FUNCTION_START(fmpz_factor_trial, state)
{
    int i, j;
    fmpz_t x;
    mpz_t y1;

    fmpz_init(x);
    mpz_init(y1);

    /* Some corner cases */
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

    /* regression test */
    {
        fmpz_factor_t fac;

        fmpz_factor_init(fac);

        fmpz_set_str(x, "-27881013806671883810", 10);

        fmpz_factor_trial(fac, x, 0);

        fmpz_factor_clear(fac);
        fmpz_factor_init(fac);

        fmpz_factor_trial(fac, x, 0);

        fmpz_factor_clear(fac);
    }

    fmpz_clear(x);

    TEST_FUNCTION_END(state);
}
#undef check
