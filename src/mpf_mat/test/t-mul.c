/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpf_mat.h"

#define MPF_MAT_MUL_BITS (40)

TEST_FUNCTION_START(mpf_mat_mul, state)
{
    mpf_mat_t A, B, C, D, E, F, G;
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong m, n, k, l;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);
        l = n_randint(state, 10);

        mpf_mat_init(A, m, n, 200);
        mpf_mat_init(B, n, k, 200);
        mpf_mat_init(C, k, l, 200);
        mpf_mat_init(D, n, l, 200);
        mpf_mat_init(E, m, k, 200);
        mpf_mat_init(F, m, l, 200);
        mpf_mat_init(G, m, l, 200);

        mpf_mat_randtest(A, state, 200);
        mpf_mat_randtest(B, state, 200);
        mpf_mat_randtest(C, state, 200);

        mpf_mat_mul(D, B, C);
        mpf_mat_mul(E, A, B);
        mpf_mat_mul(F, A, D);
        mpf_mat_mul(G, E, C);

        if (!mpf_mat_approx_equal(F, G, MPF_MAT_MUL_BITS))
        {
            flint_printf("FAIL: results not equal\n");
            mpf_mat_print(F);
            mpf_mat_print(G);
            fflush(stdout);
            flint_abort();
        }

        if (n == k)
	{
            mpf_mat_mul(A, A, B);

            if (!mpf_mat_equal(A, E))
            {
                flint_printf("FAIL: aliasing failed\n");
                fflush(stdout);
                flint_abort();
            }
        }

        mpf_mat_clear(A);
        mpf_mat_clear(B);
        mpf_mat_clear(C);
        mpf_mat_clear(D);
        mpf_mat_clear(E);
        mpf_mat_clear(F);
        mpf_mat_clear(G);
    }

    TEST_FUNCTION_END(state);
}
