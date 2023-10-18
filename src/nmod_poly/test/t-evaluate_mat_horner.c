/*
    Copyright (C) 2014 Ashish Kedia

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_evaluate_mat_horner, state)
{
    int i, j;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        nmod_mat_t A, B;
        mp_limb_t sum, n = n_randtest_not_zero(state);
        slong m, k;

        nmod_poly_init(a, n);
        nmod_poly_randtest(a, state, n_randint(state, 50));

        m = n_randint(state, 20);

        nmod_mat_init(A, m, m, n);
        nmod_mat_init(B, m, m, n);
        nmod_mat_one(A);
        nmod_poly_evaluate_mat_horner(B, a, A);

        sum = 0;
        for (j = 0; j < a->length; j++)
           sum = n_addmod(sum, nmod_poly_get_coeff_ui(a, j), n);

        for(k = 0; k < m; k++)
            nmod_mat_entry(A, k, k) = nmod_mul(nmod_mat_entry(A, k, k), sum, A->mod);

        if (!nmod_mat_equal(A, B))
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu\n", a->length, a->mod.n);
            flint_printf("sum = %wu\n", sum);
            nmod_poly_print(a), flint_printf("\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_mat_clear(A);
        nmod_mat_clear(B);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        nmod_mat_t A, B, C;
        mp_limb_t n = n_randtest_not_zero(state);
        slong m;

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, 50));
        nmod_poly_randtest(b, state, n_randint(state, 50));

        m = n_randint(state, 20);

        nmod_mat_init(A, m, m, n);
        nmod_mat_init(B, m, m, n);
        nmod_mat_init(C, m, m, n);
        nmod_mat_randtest(A, state);

        nmod_poly_evaluate_mat_horner(B, a, A);
        nmod_poly_evaluate_mat_horner(C, b, A);
        nmod_mat_add(C, B, C);

        nmod_poly_add(a, a, b);
        nmod_poly_evaluate_mat_horner(B, a, A);

        if (!nmod_mat_equal(B, C))
        {
            flint_printf("FAIL:\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
