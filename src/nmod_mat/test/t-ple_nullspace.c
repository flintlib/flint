/*
    Copyright (C) 2026 Vincent Neiger

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

TEST_FUNCTION_START(nmod_mat_ple_nullspace, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        nmod_mat_t A, At, LU, X, Xc, Y, Yt, Z;
        slong m, n, r, rank, d, i, j, k;
        slong * P, * Q, * mrp_r, * mrp_c, * rws;
        ulong p;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        if (n_randint(state, 10) == 0)
        {
            m = n_randint(state, 150);
            n = n_randint(state, 150);
        }
        p = _rp_rand_prime(state);

        nmod_mat_init(A, m, n, p);
        mrp_r = flint_malloc((FLINT_MIN(m, n) + 1) * sizeof(slong));
        mrp_c = flint_malloc((FLINT_MIN(m, n) + 1) * sizeof(slong));
        r = n_randint(state, FLINT_MIN(m, n) + 1);
        if (n_randint(state, 2))
            r = FLINT_MIN(r, 3);
        _rp_rand_positions(mrp_r, mrp_c, m, n, r, state);
        nmod_mat_randrank_profile(A, state, mrp_r, mrp_c, r);

        P = flint_malloc((m + 1) * sizeof(slong));
        Q = flint_malloc((n + 1) * sizeof(slong));
        rws = flint_malloc((m + 1) * sizeof(slong));
        nmod_mat_init_set(LU, A);
        rank = nmod_mat_ple(P, Q, LU, 0);
        if (rank != r)
            TEST_FUNCTION_FAIL("rank\n");

        /* right kernel: A X = 0, the basis from the reduced row echelon
           form, the compact form */
        d = n - rank;
        nmod_mat_init(X, n, d, p);
        nmod_mat_randtest(X, state);
        nmod_mat_ple_nullspace(X, LU, Q, rank);
        _rp_kernel_rref(Y, A);
        if (!nmod_mat_equal(X, Y))
            TEST_FUNCTION_FAIL(
                    "right kernel\nA = %{nmod_mat}\nX = %{nmod_mat}\nY = %{nmod_mat}\n",
                    A, X, Y);
        nmod_mat_init(Z, m, d, p);
        nmod_mat_mul(Z, A, X);
        if (!nmod_mat_is_zero(Z))
            TEST_FUNCTION_FAIL("A X != 0\n");
        nmod_mat_clear(Z);

        nmod_mat_init(Xc, rank, d, p);
        nmod_mat_randtest(Xc, state);
        nmod_mat_ple_nullspace_compact(Xc, LU, Q, rank);
        for (j = 0; j < rank; j++)
            for (k = 0; k < d; k++)
                if (nmod_mat_entry(Xc, j, k) != nmod_mat_entry(X, Q[j], k))
                    TEST_FUNCTION_FAIL("right kernel, compact form\n");
        nmod_mat_clear(Xc);
        nmod_mat_clear(X);
        nmod_mat_clear(Y);

        /* left kernel: X A = 0, the basis from the reduced row echelon
           form of the transpose, the compact form */
        d = m - rank;
        nmod_mat_init(X, d, m, p);
        nmod_mat_randtest(X, state);
        nmod_mat_ple_left_nullspace(X, LU, P, rank);
        nmod_mat_init(At, n, m, p);
        nmod_mat_transpose(At, A);
        _rp_kernel_rref(Yt, At);
        nmod_mat_init(Y, d, m, p);
        nmod_mat_transpose(Y, Yt);
        if (!nmod_mat_equal(X, Y))
            TEST_FUNCTION_FAIL(
                    "left kernel\nA = %{nmod_mat}\nX = %{nmod_mat}\nY = %{nmod_mat}\n",
                    A, X, Y);
        nmod_mat_init(Z, d, n, p);
        nmod_mat_mul(Z, X, A);
        if (!nmod_mat_is_zero(Z))
            TEST_FUNCTION_FAIL("X A != 0\n");
        nmod_mat_clear(Z);

        nmod_mat_init(Xc, d, rank, p);
        nmod_mat_randtest(Xc, state);
        nmod_mat_ple_left_nullspace_compact(Xc, rws, LU, P, rank);
        /* rws: the row rank profile in increasing order, then the others */
        flint_sort(mrp_r, r, sizeof(slong), _rp_slong_cmp, NULL);
        for (i = 0; i < m; i++)
            if ((i < rank && rws[i] != mrp_r[i]) || (i > rank && rws[i] <= rws[i - 1])
                    || (i >= rank && rws[i] != P[i]))
                TEST_FUNCTION_FAIL("left kernel, compact form: rows\n");
        for (k = 0; k < d; k++)
            for (j = 0; j < rank; j++)
                if (nmod_mat_entry(Xc, k, j) != nmod_mat_entry(X, k, rws[j]))
                    TEST_FUNCTION_FAIL("left kernel, compact form\n");

        nmod_mat_clear(Xc);
        nmod_mat_clear(X);
        nmod_mat_clear(Y);
        nmod_mat_clear(Yt);
        nmod_mat_clear(At);
        nmod_mat_clear(A);
        nmod_mat_clear(LU);
        flint_free(P);
        flint_free(Q);
        flint_free(rws);
        flint_free(mrp_r);
        flint_free(mrp_c);
    }

    TEST_FUNCTION_END(state);
}
