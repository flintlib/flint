/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "rank_profile_helpers.h"

TEST_FUNCTION_START(nmod_mat_nullspace, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, ker;
        ulong mod;
        slong m, n, d, r, nullity, nulrank;

        m = n_randint(state, 30);
        n = n_randint(state, 30);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            mod = n_randtest_prime(state, 0);
            d = n_randint(state, 2*m*n + 1);

            nmod_mat_init(A, m, n, mod);
            nmod_mat_init(ker, n, n, mod);
            nmod_mat_init(B, m, n, mod);

            nmod_mat_randrank(A, state, r);
            /* Densify */
            if (n_randlimb(state) % 2)
                nmod_mat_randops(A, state, d);

            nullity = nmod_mat_nullspace(ker, A);
            nulrank = nmod_mat_rank(ker);

            if (nullity != nulrank)
                TEST_FUNCTION_FAIL(
                        "rank(ker) != nullity\n"
                        "A = %{nmod_mat}\n",
                        A);

            if (nullity + r != n)
                TEST_FUNCTION_FAIL(
                        "nullity + rank != n\n"
                        "A = %{nmod_mat}\n",
                        A);

            nmod_mat_mul(B, A, ker);

            if (nmod_mat_rank(B) != 0)
                TEST_FUNCTION_FAIL(
                        "A * ker != 0\n"
                        "A = %{nmod_mat}\n",
                        A);

            /* the basis with identity at the columns without pivot (as
               computed from the reduced row echelon form), the other
               columns of ker zero */
            {
                nmod_mat_t K, W;
                slong j, k;

                _rp_kernel_rref(K, A);
                nmod_mat_window_init(W, ker, 0, 0, n, nullity);
                if (!nmod_mat_equal(W, K))
                    TEST_FUNCTION_FAIL(
                            "basis\n"
                            "A = %{nmod_mat}\n"
                            "ker = %{nmod_mat}\n"
                            "K = %{nmod_mat}\n",
                            A, ker, K);
                nmod_mat_window_clear(W);
                for (j = 0; j < n; j++)
                    for (k = nullity; k < n; k++)
                        if (nmod_mat_entry(ker, j, k) != 0)
                            TEST_FUNCTION_FAIL("nonzero entry beyond the basis\n");

                /* aliasing */
                if (m == n)
                {
                    nmod_mat_t A2;
                    nmod_mat_init_set(A2, A);
                    if (nmod_mat_nullspace(A2, A2) != nullity || !nmod_mat_equal(A2, ker))
                        TEST_FUNCTION_FAIL("aliasing\n");
                    nmod_mat_clear(A2);
                }

                nmod_mat_clear(K);
            }

            nmod_mat_clear(A);
            nmod_mat_clear(ker);
            nmod_mat_clear(B);
        }
    }

    TEST_FUNCTION_END(state);
}
