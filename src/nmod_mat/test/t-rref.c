/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "perm.h"
#include "nmod.h"
#include "nmod_mat.h"

static int check_rref_form(slong * FLINT_UNUSED(perm), nmod_mat_t A, slong rank)
{
    slong i, j, k, prev_pivot;

    /* bottom should be zero */
    for (i = rank; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (nmod_mat_entry(A, i, j) != 0)
                return 0;

    prev_pivot = -1;

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (nmod_mat_entry(A, i, j) != 0)
            {
                /* pivot should have a higher column index than previous */
                if (j <= prev_pivot)
                    return 0;

                /* column should be 0 ... 0 1 0 ... 0 */
                for (k = 0; k < rank; k++)
                    if (nmod_mat_entry(A, k, j) != (i == k))
                        return 0;

                prev_pivot = j;
                break;
            }
        }
    }

    return 1;
}

TEST_FUNCTION_START(nmod_mat_rref, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        ulong mod;
        slong j, k, m, n, rank1, rank2;
        slong *perm;
        int equal;
        ulong c;

        mod = n_randtest_prime(state, 0);

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        perm = _perm_init(2*m);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(D, 2*m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_init_set(B, A);
        nmod_mat_init_set(C, A);

        rank1 = nmod_mat_rref(B);

        if (!check_rref_form(perm, B, rank1))
            TEST_FUNCTION_FAIL(
                    "Malformed rref\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n",
                    A, B);

        /* Concatenate the original matrix with the rref, scramble the rows,
           and check that the rref is the same */
        _perm_randtest(perm, 2 * m, state);

        for (j = 0; j < m; j++)
        {
            do { c = n_randint(state, mod); } while (c == 0);
            for (k = 0; k < n; k++)
                nmod_mat_entry(D, perm[j], k) =
                    nmod_mul(nmod_mat_entry(A, j, k), c, A->mod);
        }

        for (j = 0; j < m; j++)
        {
            do { c = n_randint(state, mod); } while (c == 0);
            for (k = 0; k < n; k++)
                nmod_mat_entry(D, perm[m + j], k) =
                    nmod_mul(nmod_mat_entry(B, j, k), c, A->mod);
        }

        rank2 = nmod_mat_rref(D);
        equal = (rank1 == rank2);

        if (equal)
        {
            for (j = 0; j < rank2; j++)
                for (k = 0; k < n; k++)
                    equal = equal && (nmod_mat_entry(B, j, k) ==
                                        nmod_mat_entry(D, j, k));
            for (j = rank2; j < 2 * rank2; j++)
                for (k = 0; k < n; k++)
                    equal = equal && (nmod_mat_entry(D, j, k) == 0);
        }

        if (!equal)
            TEST_FUNCTION_FAIL(
                    "rank1 = %wd, rank2 = %wd\n"
                    "A = %{nmod_mat}\n"
                    "B = %{nmod_mat}\n"
                    "D = %{nmod_mat}\n",
                    rank1, rank2, A, B, D);

        _perm_clear(perm);
        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}
