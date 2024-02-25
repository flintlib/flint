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
    int i;
    fmpz_t x;

    fmpz_init(x);

    /* Some corner cases */
    fmpz_set_si(x, COEFF_MAX);
    check(x);
    fmpz_set_si(x, COEFF_MIN);
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

    /* Large negative integers */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_set_ui(x, 10 + n_randint(state, 10));
        fmpz_pow_ui(x, x, n_randint(state, 100));
        fmpz_neg(x, x);
        check(x);
    }

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
