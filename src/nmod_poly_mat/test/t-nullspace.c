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

TEST_FUNCTION_START(nmod_poly_mat_nullspace, state)
{
    slong i;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, N, AN;
        slong n, m, deg, rank, nullity;
        float density;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 13);
        n = n_randint(state, 13);
        deg = 1 + n_randint(state, 5);
        density = n_randint(state, 100) * 0.01;

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(N, n, n, mod);
        nmod_poly_mat_init(AN, m, n, mod);

        nmod_poly_mat_randtest_sparse(A, state, deg, density);

        rank = nmod_poly_mat_rank(A);
        nullity = nmod_poly_mat_nullspace(N, A);

        if (nullity + rank != n)
        {
            flint_printf("FAIL: wrong nullity!\n");
            flint_printf("rank = %wd\n", rank);
            flint_printf("nullity = %wd\n", nullity);
            nmod_poly_mat_print(A, "x");
            flint_printf("\n");
            nmod_poly_mat_print(N, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        if (nmod_poly_mat_rank(N) != nullity)
        {
            flint_printf("FAIL: wrong rank(N) != nullity!\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_mul(AN, A, N);

        if (!nmod_poly_mat_is_zero(AN))
        {
            flint_printf("FAIL: A * N != 0\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(N);
        nmod_poly_mat_clear(AN);
    }

    TEST_FUNCTION_END(state);
}
