/*
    Copyright (C) 2010, 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
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
void check(slong * P, nmod_mat_t LU, const nmod_mat_t A, slong rank)
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
    {
        for (j = i; j < n; j++)
        {
            if (nmod_mat_entry(LU, i, j) != 0)
            {
                flint_printf("FAIL: wrong shape!\n");
                fflush(stdout);
                flint_abort();
            }
        }
    }

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
    {
        flint_printf("FAIL\n");
        flint_printf("A:\n");
        nmod_mat_print_pretty(A);
        flint_printf("LU:\n");
        nmod_mat_print_pretty(LU);
        flint_printf("B:\n");
        nmod_mat_print_pretty(B);
        fflush(stdout);
        flint_abort();
    }

    nmod_mat_clear(B);
    nmod_mat_clear(L);
    nmod_mat_clear(U);
}
#endif

TEST_FUNCTION_START(nmod_mat_lu_classical, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, LU;
        mp_limb_t mod;
        slong m, n, r, d, rank;
        slong * P;

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
                nmod_mat_randops(A, d, state);
            }

            nmod_mat_init_set(LU, A);
            P = flint_malloc(sizeof(slong) * m);

            rank = nmod_mat_lu_classical(P, LU, 0);

            if (r != rank)
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                flint_printf("A:");
                nmod_mat_print_pretty(A);
                flint_printf("LU:");
                nmod_mat_print_pretty(LU);
                fflush(stdout);
                flint_abort();
            }

            check(P, LU, A, rank);

            nmod_mat_clear(A);
            nmod_mat_clear(LU);
            flint_free(P);
        }
    }

    TEST_FUNCTION_END(state);
}
