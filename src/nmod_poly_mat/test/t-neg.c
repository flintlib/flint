/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_neg, state)
{
    slong i;

    /* Check evaluation homomorphism */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B;
        nmod_mat_t a, b, c;
        mp_limb_t x, mod;
        slong m, n, deg;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(B, m, n, mod);

        nmod_mat_init(a, m, n, mod);
        nmod_mat_init(b, m, n, mod);
        nmod_mat_init(c, m, n, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_neg(B, A);

        x = n_randint(state, mod);

        nmod_poly_mat_evaluate_nmod(a, A, x);
        nmod_poly_mat_evaluate_nmod(b, B, x);
        nmod_mat_neg(c, a);

        if (!nmod_mat_equal(b, c))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            nmod_poly_mat_print(A, "x");
            flint_printf("B:\n");
            nmod_poly_mat_print(B, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);

        nmod_mat_clear(a);
        nmod_mat_clear(b);
        nmod_mat_clear(c);
    }

    /* Check aliasing B and A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B;
        slong m, n, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(B, m, n, mod);

        nmod_poly_mat_randtest(A, state, deg);

        nmod_poly_mat_neg(B, A);
        nmod_poly_mat_neg(A, A);

        if (!nmod_poly_mat_equal(B, A))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            nmod_poly_mat_print(A, "x");
            flint_printf("B:\n");
            nmod_poly_mat_print(B, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
