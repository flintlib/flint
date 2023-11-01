/*
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_get_set_coeff_mat, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t pmat1;
        nmod_poly_mat_t pmat2;
        nmod_mat_t cmat1;
        nmod_mat_t cmat2;
        mp_limb_t mod;
        slong m, n, deg;
        int jx;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 3 + n_randint(state, 10);

        nmod_poly_mat_init(pmat1, m, n, mod);
        nmod_poly_mat_init(pmat2, m, n, mod);
        nmod_mat_init(cmat1, m, n, mod);
        nmod_mat_init(cmat2, m, n, mod);

        // test 1: set then get does not change values
        nmod_mat_randtest(cmat1, state);
        nmod_poly_mat_set_coeff_mat(pmat1, cmat1, deg-1);
        nmod_poly_mat_get_coeff_mat(cmat2, pmat1, deg-1);
        if (!nmod_mat_equal(cmat1, cmat2))
        {
            flint_printf("FAIL (set then get):\n");
            flint_printf("pmat1:\n");
            nmod_poly_mat_print(pmat1, "x");
            flint_printf("cmat1:\n");
            nmod_mat_print(cmat1);
            flint_printf("cmat2:\n");
            nmod_mat_print(cmat2);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        // test 2: copying via repeated get-set works
        nmod_poly_mat_randtest(pmat1, state, deg);
        for (jx = 0; jx < nmod_poly_mat_max_length(pmat1); jx++)
        {
            nmod_poly_mat_get_coeff_mat(cmat1, pmat1, jx);
            nmod_poly_mat_set_coeff_mat(pmat2, cmat1, jx);
        }

        if (!nmod_poly_mat_equal(pmat1, pmat2))
        {
            flint_printf("FAIL (simulate copy):\n");
            flint_printf("pmat1:\n");
            nmod_poly_mat_print(pmat1, "x");
            flint_printf("pmat2:\n");
            nmod_poly_mat_print(pmat2, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_clear(pmat1);
        nmod_poly_mat_clear(pmat2);
        nmod_mat_clear(cmat1);
        nmod_mat_clear(cmat2);
    }

    TEST_FUNCTION_END(state);
}
