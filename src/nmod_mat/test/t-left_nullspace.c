/*
    Copyright (C) 2026 Edgar Costa

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

TEST_FUNCTION_START(nmod_mat_left_nullspace, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, ker;
        ulong mod;
        slong m, n, d, r, nullity, nulrank;

        m = n_randint(state, 30);
        n = n_randint(state, 30);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            mod = n_randtest_prime(state, 0);
            d = n_randint(state, 2 * m * n + 1);

            nmod_mat_init(A, m, n, mod);
            nmod_mat_init(ker, m, m, mod);
            nmod_mat_init(B, m, n, mod);

            nmod_mat_randrank(A, state, r);
            if (n_randlimb(state) % 2)
                nmod_mat_randops(A, state, d);

            nullity = nmod_mat_left_nullspace(ker, A);
            nulrank = nmod_mat_rank(ker);

            if (nullity != nulrank)
                TEST_FUNCTION_FAIL(
                        "rank(ker) != nullity\n"
                        "A = %{nmod_mat}\n",
                        A);

            if (nullity + r != m)
                TEST_FUNCTION_FAIL(
                        "nullity + rank != m\n"
                        "A = %{nmod_mat}\n",
                        A);

            /* Rows of ker are left nullspace vectors: ker * A = 0. */
            nmod_mat_mul(B, ker, A);

            if (nmod_mat_rank(B) != 0)
                TEST_FUNCTION_FAIL(
                        "ker * A != 0\n"
                        "A = %{nmod_mat}\n",
                        A);

            /* the basis with identity at the rows of A not in its row rank
               profile (the transpose of the basis of the right kernel of
               the transpose, from its reduced row echelon form), the
               other rows of ker zero */
            {
                nmod_mat_t At, K, Kt, W;
                slong j, k;

                nmod_mat_init(At, n, m, mod);
                nmod_mat_transpose(At, A);
                _rp_kernel_rref(Kt, At);
                nmod_mat_init(K, nullity, m, mod);
                nmod_mat_transpose(K, Kt);
                nmod_mat_window_init(W, ker, 0, 0, nullity, m);
                if (!nmod_mat_equal(W, K))
                    TEST_FUNCTION_FAIL(
                            "basis\n"
                            "A = %{nmod_mat}\n"
                            "ker = %{nmod_mat}\n"
                            "K = %{nmod_mat}\n",
                            A, ker, K);
                nmod_mat_window_clear(W);
                for (j = nullity; j < m; j++)
                    for (k = 0; k < m; k++)
                        if (nmod_mat_entry(ker, j, k) != 0)
                            TEST_FUNCTION_FAIL("nonzero entry beyond the basis\n");

                nmod_mat_clear(At);
                nmod_mat_clear(K);
                nmod_mat_clear(Kt);
            }

            nmod_mat_clear(A);
            nmod_mat_clear(ker);
            nmod_mat_clear(B);
        }
    }

    TEST_FUNCTION_END(state);
}
