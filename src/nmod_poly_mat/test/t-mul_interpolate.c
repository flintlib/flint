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

TEST_FUNCTION_START(nmod_poly_mat_mul_interpolate, state)
{
    slong i;

    /* Check evaluation homomorphism */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, C;
        nmod_mat_t a, b, c, d;
        mp_limb_t mod, x;
        slong m, n, k, deg;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        k = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(B, n, k, mod);
        nmod_poly_mat_init(C, m, k, mod);

        nmod_mat_init(a, m, n, mod);
        nmod_mat_init(b, n, k, mod);
        nmod_mat_init(c, m, k, mod);
        nmod_mat_init(d, m, k, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_randtest(B, state, deg);
        nmod_poly_mat_randtest(C, state, deg);  /* noise in output */

        if (nmod_poly_mat_max_length(A)
            + nmod_poly_mat_max_length(B) - 1 <= mod)
        {
            nmod_poly_mat_mul_interpolate(C, A, B);

            x = n_randint(state, mod);

            nmod_poly_mat_evaluate_nmod(a, A, x);
            nmod_poly_mat_evaluate_nmod(b, B, x);
            nmod_poly_mat_evaluate_nmod(d, C, x);
            nmod_mat_mul(c, a, b);

            if (!nmod_mat_equal(c, d))
            {
                flint_printf("FAIL:\n");
                flint_printf("A:\n");
                nmod_poly_mat_print(A, "x");
                flint_printf("B:\n");
                nmod_poly_mat_print(B, "x");
                flint_printf("C:\n");
                nmod_poly_mat_print(C, "x");
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);

        nmod_mat_clear(a);
        nmod_mat_clear(b);
        nmod_mat_clear(c);
        nmod_mat_clear(d);
    }

    /* Check aliasing C and A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, C;
        slong m, n, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(B, n, n, mod);
        nmod_poly_mat_init(C, m, n, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_randtest(B, state, deg);
        nmod_poly_mat_randtest(C, state, deg);  /* noise in output */

        if (nmod_poly_mat_max_length(A)
            + nmod_poly_mat_max_length(B) - 1 <= mod)
        {

            nmod_poly_mat_mul_interpolate(C, A, B);
            nmod_poly_mat_mul_interpolate(A, A, B);

            if (!nmod_poly_mat_equal(C, A))
            {
                flint_printf("FAIL:\n");
                flint_printf("A:\n");
                nmod_poly_mat_print(A, "x");
                flint_printf("B:\n");
                nmod_poly_mat_print(B, "x");
                flint_printf("C:\n");
                nmod_poly_mat_print(C, "x");
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
    }

    /* Check aliasing C and B */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, C;
        slong m, n, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(A, m, m, mod);
        nmod_poly_mat_init(B, m, n, mod);
        nmod_poly_mat_init(C, m, n, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_randtest(B, state, deg);
        nmod_poly_mat_randtest(C, state, deg);  /* noise in output */

        if (nmod_poly_mat_max_length(A)
            + nmod_poly_mat_max_length(B) - 1 <= mod)
        {
            nmod_poly_mat_mul_interpolate(C, A, B);
            nmod_poly_mat_mul_interpolate(B, A, B);

            if (!nmod_poly_mat_equal(C, B))
            {
                flint_printf("FAIL:\n");
                flint_printf("A:\n");
                nmod_poly_mat_print(A, "x");
                flint_printf("B:\n");
                nmod_poly_mat_print(B, "x");
                flint_printf("C:\n");
                nmod_poly_mat_print(C, "x");
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
