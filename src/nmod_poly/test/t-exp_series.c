/*
    Copyright (C) 2010 William Hart
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

TEST_FUNCTION_START(nmod_poly_exp_series, state)
{
    int i, result = 1;

    /* Check exp(A+B) = exp(A) * exp(B) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, B, AB, expA, expB, expAB, S;
        slong n;
        slong N = 100;
        mp_limb_t mod;

        /* Make sure to workout the Newton code */
        if (n_randint(state, 100) == 1)
            N = 20000;

        mod = n_randtest_prime(state, 0);
        n = n_randtest(state) % N;
        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(B, mod);
        nmod_poly_init(AB, mod);
        nmod_poly_init(expA, mod);
        nmod_poly_init(expB, mod);
        nmod_poly_init(expAB, mod);
        nmod_poly_init(S, mod);

        nmod_poly_randtest(A, state, n_randint(state, N));
        nmod_poly_set_coeff_ui(A, 0, UWORD(0));
        nmod_poly_randtest(B, state, n_randint(state, N));
        nmod_poly_set_coeff_ui(B, 0, UWORD(0));

        /* Randomly generate a monomial */
        if (n_randlimb(state) % 100 == 0)
        {
            nmod_poly_zero(A);
            nmod_poly_set_coeff_ui(A, n_randlimb(state) % (n+5), \
                n_randtest_not_zero(state) % mod);
            nmod_poly_set_coeff_ui(A, 0, UWORD(0));
        }

        nmod_poly_exp_series(expA, A, n);
        nmod_poly_exp_series(expB, B, n);
        nmod_poly_add(AB, A, B);
        nmod_poly_exp_series(expAB, AB, n);
        nmod_poly_mullow(S, expA, expB, n);

        result = nmod_poly_equal(S, expAB);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd, mod = %wu\n", n, mod);
            flint_printf("A: "); nmod_poly_print(A), flint_printf("\n\n");
            flint_printf("B: "); nmod_poly_print(B), flint_printf("\n\n");
            flint_printf("exp(A): "); nmod_poly_print(expA), flint_printf("\n\n");
            flint_printf("exp(B): "); nmod_poly_print(expB), flint_printf("\n\n");
            flint_printf("exp(A+B):       "); nmod_poly_print(expAB), flint_printf("\n\n");
            flint_printf("exp(A)*exp(B): "); nmod_poly_print(S), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(A);
        nmod_poly_clear(B);
        nmod_poly_clear(AB);
        nmod_poly_clear(expA);
        nmod_poly_clear(expB);
        nmod_poly_clear(expAB);
        nmod_poly_clear(S);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, B;
        slong n;
        mp_limb_t mod;
        mod = n_randtest_prime(state, 0);

        if (i < 30)
            n = n_randint(state, 20000);
        else
            n = n_randtest(state) % 50;

        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(B, mod);
        nmod_poly_randtest(A, state, n);
        nmod_poly_set_coeff_ui(A, 0, UWORD(0));

        nmod_poly_exp_series(B, A, n);
        nmod_poly_exp_series(A, A, n);

        result = nmod_poly_equal(A, B);
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(A), flint_printf("\n\n");
            nmod_poly_print(B), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(A);
        nmod_poly_clear(B);
    }

    TEST_FUNCTION_END(state);
}
