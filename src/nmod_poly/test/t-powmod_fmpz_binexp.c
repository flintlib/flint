/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "fmpz.h"

TEST_FUNCTION_START(nmod_poly_powmod_fmpz_binexp, state)
{
    int i, result;

    /* Aliasing of res and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, t, f;
        mp_limb_t n;
        fmpz_t exp;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);;

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_powmod_fmpz_binexp(res1, a, exp, f);
        nmod_poly_powmod_fmpz_binexp(a, a, exp, f);

        result = (nmod_poly_equal(res1, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp);
        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    /* Aliasing of res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, t, f;
        mp_limb_t n;
        fmpz_t exp;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_powmod_fmpz_binexp(res1, a, exp, f);
        nmod_poly_powmod_fmpz_binexp(f, a, exp, f);

        result = (nmod_poly_equal(res1, f));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp);
        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(res1);
        nmod_poly_clear(t);
    }

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, res1, res2, t, f;
        mp_limb_t n;
        fmpz_t exp;
        int j;

        fmpz_init(exp);

        n = n_randtest_prime(state, 0);
        fmpz_randtest_unsigned(exp, state, n_randint(state, 100) + 1);

        nmod_poly_init(a, n);
        nmod_poly_init(f, n);
        nmod_poly_init(res1, n);
        nmod_poly_init(res2, n);
        nmod_poly_init(t, n);

        nmod_poly_randtest(a, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(f, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(f));

        nmod_poly_powmod_fmpz_binexp(res1, a, exp, f);

	if (fmpz_cmp_ui(exp, 32) <= 0)
        {
           nmod_poly_zero(res2);
           if (nmod_poly_length(f) > 1)
	       nmod_poly_set_coeff_ui(res2, 0, 1);
           for (j = 1; j <= fmpz_get_ui(exp); j++)
               nmod_poly_mulmod(res2, res2, a, f);
	} else
        {
	   fmpz_sub_ui(exp, exp, 1);
           nmod_poly_powmod_fmpz_binexp(res2, a, exp, f);
           nmod_poly_mulmod(res2, res2, a, f);
	   fmpz_add_ui(exp, exp, 1);
	}

        result = (nmod_poly_equal(res1, res2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: "); fmpz_print(exp); flint_printf("\n\n");
            flint_printf("a:\n"); nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("res1:\n"); nmod_poly_print(res1), flint_printf("\n\n");
            flint_printf("res2:\n"); nmod_poly_print(res2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(exp);
        nmod_poly_clear(a);
        nmod_poly_clear(f);
        nmod_poly_clear(res1);
        nmod_poly_clear(res2);
        nmod_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
