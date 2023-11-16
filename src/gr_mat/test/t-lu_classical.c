/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "gr_mat.h"

/* Defined in t-lu.c, t-lu_classical.c, t-lu_recursive.c */
#ifndef perm
#define perm perm
void perm(gr_mat_t A, slong * P)
{
    slong i;
    gr_ptr * tmp;

    if (A->c == 0 || A->r == 0)
        return;

    tmp = flint_malloc(sizeof(gr_ptr) * A->r);

    for (i = 0; i < A->r; i++) tmp[P[i]] = A->rows[i];
    for (i = 0; i < A->r; i++) A->rows[i] = tmp[i];

    flint_free(tmp);
}
#endif

/* Defined in t-lu.c, t-lu_classical.c, t-lu_recursive.c */
#ifndef check
#define check check
void check(slong * P, gr_mat_t LU, const gr_mat_t A, slong rank, gr_ctx_t ctx)
{
    gr_mat_t B, L, U;
    slong m, n, i, j;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    m = A->r;
    n = A->c;

    gr_mat_init(B, m, n, ctx);
    gr_mat_init(L, m, m, ctx);
    gr_mat_init(U, m, n, ctx);

    rank = FLINT_ABS(rank);

    for (i = rank; i < FLINT_MIN(m, n); i++)
    {
        for (j = i; j < n; j++)
        {
            if (gr_is_zero(GR_MAT_ENTRY(LU, i, j, sz), ctx) != T_TRUE)
            {
                flint_printf("FAIL: wrong shape!\n");
                flint_abort();
            }
        }
    }

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < FLINT_MIN(i, n); j++)
            status |= gr_set(GR_MAT_ENTRY(L, i, j, sz), GR_MAT_ENTRY(LU, i, j, sz), ctx);
        if (i < rank)
            status |= gr_one(GR_MAT_ENTRY(L, i, i, sz), ctx);
        for (j = i; j < n; j++)
            status |= gr_set(GR_MAT_ENTRY(U, i, j, sz), GR_MAT_ENTRY(LU, i, j, sz), ctx);
    }

    status |= gr_mat_mul(B, L, U, ctx);
    perm(B, P);

    if (status == GR_SUCCESS && gr_mat_equal(A, B, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n");
        flint_printf("A:\n");
        gr_mat_print(A, ctx);
        flint_printf("LU:\n");
        gr_mat_print(LU, ctx);
        flint_printf("B:\n");
        gr_mat_print(B, ctx);
        flint_abort();
    }

    gr_mat_clear(B, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
}
#endif

/* Defined in t-lu.c, t-lu_classical.c, t-lu_recursive.c */
#ifndef _gr_mat_randrank
#define _gr_mat_randrank _gr_mat_randrank
void
_gr_mat_randrank(gr_mat_t mat, flint_rand_t state, slong rank, slong bits, gr_ctx_t ctx)
{
    fmpz_mat_t A;
    fmpz_mat_init(A, mat->r, mat->c);
    fmpz_mat_randrank(A, state, rank, bits);
    GR_MUST_SUCCEED(gr_mat_set_fmpz_mat(mat, A, ctx));
    fmpz_mat_clear(A);
}
#endif

TEST_FUNCTION_START(gr_mat_lu_classical, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, LU;
        slong m, n, r, d, rank;
        slong * P;
        int status;

        if (n_randint(state, 2))
            gr_ctx_init_fmpz(ctx);
        else
            gr_ctx_init_fmpq(ctx);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            gr_mat_init(A, m, n, ctx);
            gr_mat_init(LU, m, n, ctx);

            _gr_mat_randrank(A, state, r, 5, ctx);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*m*n + 1);
                GR_MUST_SUCCEED(gr_mat_randops(A, state, d, ctx));
            }

            P = flint_malloc(sizeof(slong) * m);

            status = gr_mat_lu_classical(&rank, P, LU, A, 0, ctx);

            if (status == GR_SUCCESS)
            {
                if (r != rank)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("wrong rank!\n");
                    flint_printf("A:");
                    gr_mat_print(A, ctx);
                    flint_printf("LU:");
                    gr_mat_print(LU, ctx);
                    flint_abort();
                }

                check(P, LU, A, rank, ctx);
            }

            gr_mat_clear(A, ctx);
            gr_mat_clear(LU, ctx);
            flint_free(P);
        }

        gr_ctx_clear(ctx);
    }

    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, LU;
        slong m, n, r, d, rank;
        slong * P;
        int status;

        gr_ctx_init_random(ctx, state);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            gr_mat_init(A, m, n, ctx);
            gr_mat_init(LU, m, n, ctx);

            _gr_mat_randrank(A, state, r, 5, ctx);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*m*n + 1);
                GR_MUST_SUCCEED(gr_mat_randops(A, state, d, ctx));
            }

            P = flint_malloc(sizeof(slong) * m);

            status = gr_mat_lu_classical(&rank, P, LU, A, 0, ctx);

            if (status == GR_SUCCESS)
            {
                if (rank > r)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("wrong rank!\n");
                    flint_printf("A:");
                    gr_mat_print(A, ctx);
                    flint_printf("LU:");
                    gr_mat_print(LU, ctx);
                    flint_abort();
                }

                check(P, LU, A, rank, ctx);
            }

            gr_mat_clear(A, ctx);
            gr_mat_clear(LU, ctx);
            flint_free(P);
        }

        gr_ctx_clear(ctx);
    }

    /* Test rank check for square matrices */
    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, LU;
        slong n, r, d, rank;
        slong * P;
        int status;

        if (n_randint(state, 2))
            gr_ctx_init_fmpz(ctx);
        else
            gr_ctx_init_fmpq(ctx);

        n = n_randint(state, 10);

        for (r = 0; r <= n; r++)
        {
            gr_mat_init(A, n, n, ctx);
            gr_mat_init(LU, n, n, ctx);

            _gr_mat_randrank(A, state, r, 5, ctx);

            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*n*n + 1);
                GR_MUST_SUCCEED(gr_mat_randops(A, state, d, ctx));
            }

            P = flint_malloc(sizeof(slong) * n);

            status = gr_mat_lu_classical(&rank, P, LU, A, 1, ctx);

            if (status == GR_SUCCESS)
            {
                if ((rank == 0) != (r < n || n == 0))
                {
                    flint_printf("FAIL:\n");
                    gr_ctx_println(ctx);
                    flint_printf("rank check\n");
                    flint_printf("r = %wd, rank = %wd\n", r, rank);
                    flint_printf("A:\n");
                    gr_mat_print(A, ctx);
                    flint_printf("LU:\n");
                    gr_mat_print(LU, ctx);
                    flint_abort();
                }

                if (rank != 0)
                    check(P, LU, A, rank, ctx);
            }

            gr_mat_clear(A, ctx);
            gr_mat_clear(LU, ctx);
            flint_free(P);
        }

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
