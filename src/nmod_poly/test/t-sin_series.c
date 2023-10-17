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

TEST_FUNCTION_START(nmod_poly_sin_series, state)
{
    int i, result = 1;

    /* Check asin(sin(A)) = A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, sinA, B;
        slong n;
        mp_limb_t mod;

        do { mod = n_randtest_prime(state, 0); } while (mod == 2);
        n = 1 + n_randtest(state) % 100;
        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(sinA, mod);
        nmod_poly_init(B, mod);

        nmod_poly_randtest(A, state, n_randint(state, 100));
        nmod_poly_set_coeff_ui(A, 0, UWORD(0));

        nmod_poly_sin_series(sinA, A, n);
        nmod_poly_asin_series(B, sinA, n);

        nmod_poly_truncate(A, n);

        result = nmod_poly_equal(A, B);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd, mod = %wu\n", n, mod);
            flint_printf("A: "); nmod_poly_print(A), flint_printf("\n\n");
            flint_printf("sin(A): "); nmod_poly_print(sinA), flint_printf("\n\n");
            flint_printf("B: "); nmod_poly_print(B), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(A);
        nmod_poly_clear(sinA);
        nmod_poly_clear(B);
    }

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t A, B;
        slong n;
        mp_limb_t mod;
        do { mod = n_randtest_prime(state, 0); } while (mod == 2);
        n = n_randtest(state) % 50;
        n = FLINT_MIN(n, mod);

        nmod_poly_init(A, mod);
        nmod_poly_init(B, mod);
        nmod_poly_randtest(A, state, n_randint(state, 50));
        nmod_poly_set_coeff_ui(A, 0, UWORD(0));

        nmod_poly_sin_series(B, A, n);
        nmod_poly_sin_series(A, A, n);

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
