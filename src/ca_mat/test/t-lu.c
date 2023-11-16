/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "ca_mat.h"

/* Defined in t-lu.c, t-lu_classical.c and t-lu_recursive.c */
#ifndef perm
#define perm perm
void perm(ca_mat_t A, slong * P)
{
    slong i;
    ca_struct ** tmp;

    if (A->c == 0 || A->r == 0)
        return;

    tmp = flint_malloc(sizeof(ca_struct *) * A->r);

    for (i = 0; i < A->r; i++) tmp[P[i]] = A->rows[i];
    for (i = 0; i < A->r; i++) A->rows[i] = tmp[i];

    flint_free(tmp);
}
#endif

/* Defined in t-lu.c, t-lu_classical.c and t-lu_recursive.c */
#ifndef check
#define check check
void check(slong * P, ca_mat_t LU, const ca_mat_t A, slong rank, ca_ctx_t ctx)
{
    ca_mat_t B, L, U;
    slong m, n, i, j;

    m = A->r;
    n = A->c;

    ca_mat_init(B, m, n, ctx);
    ca_mat_init(L, m, m, ctx);
    ca_mat_init(U, m, n, ctx);

    rank = FLINT_ABS(rank);

    for (i = rank; i < FLINT_MIN(m, n); i++)
    {
        for (j = i; j < n; j++)
        {
            if (ca_check_is_zero(ca_mat_entry(LU, i, j), ctx) != T_TRUE)
            {
                flint_printf("FAIL: wrong shape!\n");
                flint_abort();
            }
        }
    }

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < FLINT_MIN(i, n); j++)
            ca_set(ca_mat_entry(L, i, j), ca_mat_entry(LU, i, j), ctx);
        if (i < rank)
            ca_one(ca_mat_entry(L, i, i), ctx);
        for (j = i; j < n; j++)
            ca_set(ca_mat_entry(U, i, j), ca_mat_entry(LU, i, j), ctx);
    }

    ca_mat_mul(B, L, U, ctx);
    perm(B, P);

    if (ca_mat_check_equal(A, B, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n");
        flint_printf("A:\n");
        ca_mat_print(A, ctx);
        flint_printf("LU:\n");
        ca_mat_print(LU, ctx);
        flint_printf("B:\n");
        ca_mat_print(B, ctx);
        flint_abort();
    }

    ca_mat_clear(B, ctx);
    ca_mat_clear(L, ctx);
    ca_mat_clear(U, ctx);
}
#endif

/* Defined in t-lu.c, t-lu_classical.c and t-lu_recursive.c */
#ifndef ca_mat_randrank
#define ca_mat_randrank ca_mat_randrank
void
ca_mat_randrank(ca_mat_t mat, flint_rand_t state, slong rank, slong bits, ca_ctx_t ctx)
{
    fmpz_mat_t A;
    fmpz_mat_init(A, mat->r, mat->c);
    fmpz_mat_randrank(A, state, rank, bits);
    ca_mat_set_fmpz_mat(mat, A, ctx);
    fmpz_mat_clear(A);
}
#endif

TEST_FUNCTION_START(ca_mat_lu, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, LU;
        slong m, n, r, d, rank;
        slong * P;
        int success;

        ca_ctx_init(ctx);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            ca_mat_init(A, m, n, ctx);
            ca_mat_init(LU, m, n, ctx);

            ca_mat_randrank(A, state, r, 5, ctx);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*m*n + 1);
                ca_mat_randops(A, state, d, ctx);
            }

            P = flint_malloc(sizeof(slong) * m);

            success = ca_mat_lu(&rank, P, LU, A, 0, ctx);

            if (success)
            {
                if (r != rank)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("wrong rank!\n");
                    flint_printf("A:");
                    ca_mat_print(A, ctx);
                    flint_printf("LU:");
                    ca_mat_print(LU, ctx);
                    flint_abort();
                }

                check(P, LU, A, rank, ctx);
            }

            ca_mat_clear(A, ctx);
            ca_mat_clear(LU, ctx);
            flint_free(P);
        }

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
