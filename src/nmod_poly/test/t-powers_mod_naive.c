/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_powers_mod_naive, state)
{
    int i, result;

    /* Compare with powmod */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t f, g, pow;
        nmod_poly_struct * res;
        mp_limb_t n;
        ulong exp;
        slong j;

        n = n_randtest_prime(state, 0);
        exp = n_randint(state, 32);

        nmod_poly_init(f, n);
        nmod_poly_init(g, n);
        nmod_poly_init(pow, n);

        res = (nmod_poly_struct *) flint_malloc(exp*sizeof(nmod_poly_struct));

        for (j = 0; j < exp; j++)
            nmod_poly_init(res + j, n);

        nmod_poly_randtest(f, state, n_randint(state, 50));
        do {
            nmod_poly_randtest(g, state, n_randint(state, 50));
        } while (nmod_poly_is_zero(g));

        nmod_poly_powers_mod_naive(res, f, exp, g);

        result = 1;
        j = 0;

        if (exp > 0)
        {
            nmod_poly_one(pow);
            result = nmod_poly_equal(res + 0, pow);
        }

        for (j = 1 ; j < exp && result; j++)
        {
            nmod_poly_mulmod(pow, pow, f, g);
            result &= nmod_poly_equal(res + j, pow);
        }

        j--;

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp: %wu\n\n", exp);
            flint_printf("f:\n"); nmod_poly_print(f), flint_printf("\n\n");
            flint_printf("g:\n"); nmod_poly_print(g), flint_printf("\n\n");
            flint_printf("j: %w\n", j);
            flint_printf("pow:\n"); nmod_poly_print(pow), flint_printf("\n\n");
            flint_printf("res[%w]:\n", j); nmod_poly_print(res + j), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(pow);

        for (j = 0; j < exp; j++)
            nmod_poly_clear(res + j);

        flint_free(res);
    }

    TEST_FUNCTION_END(state);
}
