/*
    Copyright (C) 2022, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you grn redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"
#include "fmpz_mat.h"

static void perm(gr_mat_t A, slong * P)
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

static void check(slong * P, gr_mat_t LU, const gr_mat_t A, slong rank, gr_ctx_t ctx)
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
                flint_printf("FAIL: wrong shape\n");
                flint_printf("rank = %wd\n", rank);
                flint_printf("\n\nA:\n");
                gr_mat_print(A, ctx);
                flint_printf("\n\nLU:\n");
                gr_mat_print(LU, ctx);
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
        flint_printf("FAIL: LU does not match A\n");
        flint_printf("\n\nA:\n");
        gr_mat_print(A, ctx);
        flint_printf("\n\nLU:\n");
        gr_mat_print(LU, ctx);
        flint_printf("\n\nB:\n");
        gr_mat_print(B, ctx);
        flint_abort();
    }

    gr_mat_clear(B, ctx);
    gr_mat_clear(L, ctx);
    gr_mat_clear(U, ctx);
}

static void
_gr_mat_randrank(gr_mat_t mat, flint_rand_t state, slong rank, slong bits, gr_ctx_t ctx)
{
    fmpz_mat_t A;
    fmpz_mat_init(A, mat->r, mat->c);
    fmpz_mat_randrank(A, state, rank, bits);
    GR_MUST_SUCCEED(gr_mat_set_fmpz_mat(mat, A, ctx));
    fmpz_mat_clear(A);
}

FLINT_DLL extern gr_static_method_table _ca_methods;

void gr_mat_test_lu(gr_method_mat_lu_op lu_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t ctx2;
        gr_ctx_ptr ctxptr;
        gr_mat_t A, LU, LU2;
        slong m, n, d, r, rank, rank_lower_bound, rank_upper_bound;
        slong * P, * P2;
        int status;

        if (ctx == NULL)
        {
            gr_ctx_init_random(ctx2, state);
            ctxptr = ctx2;
        }
        else
            ctxptr = ctx;

        /* Hack: ca can have too much blowup */
        if (((gr_ctx_struct *) ctxptr)->methods == _ca_methods)
        {
            m = n_randint(state, FLINT_MIN(maxn, 4) + 1);
            n = n_randint(state, FLINT_MIN(maxn, 4) + 1);
        }
        else
        {
            m = n_randint(state, maxn + 1);
            n = n_randint(state, maxn + 1);
        }

        gr_mat_init(A, m, n, ctxptr);
        gr_mat_init(LU, m, n, ctxptr);
        gr_mat_init(LU2, m, n, ctxptr);

        P = flint_malloc(sizeof(slong) * m);
        P2 = flint_malloc(sizeof(slong) * m);

        /* Generate matrix with known bounds on rank */
        if (n_randint(state, 2))
        {
            GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctxptr));

            if (lu_impl != ((gr_method_mat_lu_op) gr_mat_lu_classical) &&
                GR_SUCCESS == gr_mat_lu_classical(&rank_lower_bound, P2, LU2, A, 0, ctxptr))
            {
                rank_upper_bound = rank_lower_bound;
            }
            else
            {
                rank_lower_bound = 0;
                rank_upper_bound = FLINT_MIN(m, n);
            }
        }
        else
        {
            r = n_randint(state, FLINT_MIN(m, n) + 1);
            _gr_mat_randrank(A, state, r, 5, ctxptr);
            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*m*n + 1);
                GR_MUST_SUCCEED(gr_mat_randops(A, state, d, ctxptr));
            }

            if (gr_ctx_is_finite_characteristic(ctxptr) == T_FALSE)
                rank_lower_bound = rank_upper_bound = r;
            else
            {
                rank_lower_bound = 0;
                rank_upper_bound = r;
            }
        }

        status = lu_impl(&rank, P, LU, A, 0, ctxptr);

        if (status == GR_SUCCESS)
        {
            /* Check shape of solution */
            check(P, LU, A, rank, ctxptr);

            /* Check rank */
            if (rank < rank_lower_bound || rank > rank_upper_bound)
            {
                flint_printf("FAIL\n");
                gr_ctx_println(ctxptr);
                flint_printf("wrong rank!\n");
                flint_printf("rank = %wd\n", rank);
                flint_printf("rank bounds = [%wd, %wd]\n", rank_lower_bound, rank_upper_bound);
                flint_printf("\n\nA:");
                gr_mat_print(A, ctxptr);
                flint_printf("\n\nLU:");
                gr_mat_print(LU, ctxptr);
                flint_abort();
            }
        }

        /* ----------------------------------------------------- */
        /* Rank check for square matrices                        */
        /* ----------------------------------------------------- */
        if (m == n)
        {
            status = lu_impl(&rank, P, LU, A, 1, ctxptr);

            if (status == GR_SUCCESS)
            {
                int ok;

                if (rank == 0)
                    ok = (n == 0) || (rank_lower_bound < n);
                else
                    ok = (rank_upper_bound == n);

                if (!ok)
                {
                    flint_printf("FAIL\n");
                    gr_ctx_println(ctxptr);
                    flint_printf("rank check\n");
                    flint_printf("rank = %wd\n", rank);
                    flint_printf("rank bounds = [%wd, %wd]\n", rank_lower_bound, rank_upper_bound);
                    flint_printf("\n\nA:\n");
                    gr_mat_print(A, ctxptr);
                    flint_printf("\n\nLU:\n");
                    gr_mat_print(LU, ctxptr);
                    flint_abort();
                }

                if (rank != 0)
                    check(P, LU, A, rank, ctxptr);
            }
        }

        gr_mat_clear(A, ctxptr);
        gr_mat_clear(LU, ctxptr);
        gr_mat_clear(LU2, ctxptr);
        flint_free(P);
        flint_free(P2);

        if (ctx == NULL)
            gr_ctx_clear(ctxptr);
    }
}
