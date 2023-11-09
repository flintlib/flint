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

TEST_FUNCTION_START(nmod_poly_mat_pow, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, C;
        slong m, j, exp, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 6);
        deg = 1 + n_randint(state, 6);
        exp = n_randint(state, 20);

        nmod_poly_mat_init(A, m, m, mod);
        nmod_poly_mat_init(B, m, m, mod);
        nmod_poly_mat_init(C, m, m, mod);

        nmod_poly_mat_randtest(A, state, deg);

        nmod_poly_mat_pow(B, A, exp);

        nmod_poly_mat_one(C);
        for (j = 0; j < exp; j++)
            nmod_poly_mat_mul(C, C, A);

        if (!nmod_poly_mat_equal(C, B))
        {
            flint_printf("FAIL:\n");
            flint_printf("exp = %wd\n", exp);
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

    /* Check aliasing B and A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B;
        slong m, exp, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 6);
        deg = 1 + n_randint(state, 6);
        exp = n_randint(state, 20);

        nmod_poly_mat_init(A, m, m, mod);
        nmod_poly_mat_init(B, m, m, mod);

        nmod_poly_mat_randtest(A, state, deg);

        nmod_poly_mat_pow(B, A, exp);
        nmod_poly_mat_pow(A, A, exp);

        if (!nmod_poly_mat_equal(A, B))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("exp = %wd\n", exp);
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
