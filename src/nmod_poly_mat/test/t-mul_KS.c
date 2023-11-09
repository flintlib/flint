/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_mul_KS, state)
{
    slong i;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, C, D;
        slong m, n, k, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 15);
        n = n_randint(state, 15);
        k = n_randint(state, 15);
        deg = 1 + n_randint(state, 15);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(B, n, k, mod);
        nmod_poly_mat_init(C, m, k, mod);
        nmod_poly_mat_init(D, m, k, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_randtest(B, state, deg);
        nmod_poly_mat_randtest(C, state, deg);  /* noise in output */

        nmod_poly_mat_mul_classical(C, A, B);
        nmod_poly_mat_mul_KS(D, A, B);

        if (!nmod_poly_mat_equal(C, D))
        {
            flint_printf("FAIL:\n");
            flint_printf("products don't agree!\n");
            flint_printf("A:\n");
            nmod_poly_mat_print(A, "x");
            flint_printf("B:\n");
            nmod_poly_mat_print(B, "x");
            flint_printf("C:\n");
            nmod_poly_mat_print(C, "x");
            flint_printf("D:\n");
            nmod_poly_mat_print(D, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
        nmod_poly_mat_clear(D);
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

        nmod_poly_mat_mul_KS(C, A, B);
        nmod_poly_mat_mul_KS(A, A, B);

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

        nmod_poly_mat_mul_KS(C, A, B);
        nmod_poly_mat_mul_KS(B, A, B);

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

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
