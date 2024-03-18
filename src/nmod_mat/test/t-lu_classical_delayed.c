/*
    Copyright (C) 2010, 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "perm.h"
#include "nmod_mat.h"

/* Defined in t-lu_classical.c, t-lu_classical_delayed.c and t-lu_recursive.c */
#ifndef perm
#define perm perm
void perm(nmod_mat_t A, slong * P)
{
    slong i;
    mp_ptr * tmp;

    if (A->c == 0 || A->r == 0)
        return;

    tmp = flint_malloc(sizeof(mp_ptr) * A->r);

    for (i = 0; i < A->r; i++) tmp[P[i]] = A->rows[i];
    for (i = 0; i < A->r; i++) A->rows[i] = tmp[i];

    flint_free(tmp);
}
#endif

/* Defined in t-lu_classical.c, t-lu_classical_delayed.c and t-lu_recursive.c */
#ifndef check
#define check check
int check(slong * P, nmod_mat_t LU, const nmod_mat_t A, slong rank)
{
    nmod_mat_t B, L, U;
    slong m, n, i, j;

    m = A->r;
    n = A->c;

    nmod_mat_init(B, m, n, A->mod.n);
    nmod_mat_init(L, m, m, A->mod.n);
    nmod_mat_init(U, m, n, A->mod.n);

    rank = FLINT_ABS(rank);

    for (i = rank; i < FLINT_MIN(m, n); i++)
        for (j = i; j < n; j++)
            if (nmod_mat_entry(LU, i, j) != 0)
                return 1;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < FLINT_MIN(i, n); j++)
            nmod_mat_entry(L, i, j) = nmod_mat_entry(LU, i, j);
        if (i < rank)
            nmod_mat_entry(L, i, i) = UWORD(1);
        for (j = i; j < n; j++)
            nmod_mat_entry(U, i, j) = nmod_mat_entry(LU, i, j);
    }

    nmod_mat_mul(B, L, U);
    perm(B, P);

    if (!nmod_mat_equal(A, B))
        return 2;

    nmod_mat_clear(B);
    nmod_mat_clear(L);
    nmod_mat_clear(U);

    return 0;
}
#endif

TEST_FUNCTION_START(nmod_mat_lu_classical_delayed, state)
{
    slong i;
    int result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, LU, LU2;
        mp_limb_t mod;
        slong m, n, r, d, rank, rank2;
        slong *P, *P2;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        mod = n_randtest_prime(state, 0);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            nmod_mat_init(A, m, n, mod);
            nmod_mat_randrank(A, state, r);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*m*n + 1);
                nmod_mat_randops(A, state, d);
            }

            nmod_mat_init_set(LU, A);
            P = flint_malloc(sizeof(slong) * m);
            rank = nmod_mat_lu_classical_delayed(P, LU, 0);

            if (r != rank)
                TEST_FUNCTION_FAIL(
                        "Wrong rank\n"
                        "A = %{nmod_mat}\n"
                        "LU = %{nmod_mat}\n",
                        A, LU);

            result = check(P, LU, A, rank);
            if (result != 0)
            {
                if (result == 1)
                    TEST_FUNCTION_FAIL("Wrong shape\n");
                else if (result == 2)
                    TEST_FUNCTION_FAIL(
                        "A = %{nmod_mat}\n"
                        "LU = %{nmod_mat}\n",
                        A, LU);
            }

            nmod_mat_init_set(LU2, A);
            P2 = flint_malloc(sizeof(slong) * m);
            rank2 = nmod_mat_lu_classical(P2, LU2, 0);

            if (r != rank || !nmod_mat_equal(LU, LU2) || !_perm_equal(P, P2, m))
                TEST_FUNCTION_FAIL(
                        "A = %{nmod_mat}\n"
                        "LU = %{nmod_mat}\n"
                        "LU2 = %{nmod_mat}\n",
                        A, LU, LU2);

            nmod_mat_set(LU, A);
            rank = nmod_mat_lu_classical_delayed(P, LU, 1);
            nmod_mat_set(LU2, A);
            rank2 = nmod_mat_lu_classical(P2, LU2, 1);

            if (rank != rank2 || (rank == r && (!nmod_mat_equal(LU, LU2) || !_perm_equal(P, P2, m))))
                TEST_FUNCTION_FAIL(
                        "A = %{nmod_mat}\n"
                        "r = %wd, rank = %wd, rank2 = %wd\n"
                        "LU = %{nmod_mat}\n"
                        "LU2 = %{nmod_mat}\n",
                        A, r, rank, rank2, LU, LU2);

            nmod_mat_clear(A);
            nmod_mat_clear(LU);
            flint_free(P);
            nmod_mat_clear(LU2);
            flint_free(P2);
        }
    }

    TEST_FUNCTION_END(state);
}
