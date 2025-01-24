/*
    Copyright (C) 2022, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you grn redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "fmpz_mat.h"

static void perm(gr_mat_t A, slong * P, gr_ctx_t ctx)
{
    slong i;
    gr_mat_t tmp;

    if (A->c == 0 || A->r == 0)
        return;

    gr_mat_init(tmp, A->r, A->c, ctx);

    for (i = 0; i < A->r; i++)
        GR_MUST_SUCCEED(_gr_vec_set(gr_mat_entry_ptr(tmp, P[i], 0, ctx),
                    gr_mat_entry_srcptr(A, i, 0, ctx), A->c, ctx));

    for (i = 0; i < A->r; i++)
        GR_MUST_SUCCEED(_gr_vec_set(gr_mat_entry_ptr(A, i, 0, ctx),
                    gr_mat_entry_srcptr(tmp, i, 0, ctx), A->c, ctx));

    gr_mat_clear(tmp, ctx);
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
    perm(B, P, ctx);

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
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t my_ctx;
        gr_ctx_ptr ctx;
        gr_mat_t A, LU, LU2;
        slong m, n, d, r, rank, rank_lower_bound, rank_upper_bound;
        slong * P, * P2;
        int status;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        /* Hack: ca can have too much blowup */
        if (((gr_ctx_struct *) ctx)->methods == _ca_methods)
        {
            m = n_randint(state, FLINT_MIN(maxn, 4) + 1);
            n = n_randint(state, FLINT_MIN(maxn, 4) + 1);
        }
        else
        {
            m = n_randint(state, maxn + 1);
            n = n_randint(state, maxn + 1);
        }

        gr_mat_init(A, m, n, ctx);
        gr_mat_init(LU, m, n, ctx);
        gr_mat_init(LU2, m, n, ctx);

        P = flint_malloc(sizeof(slong) * m);
        P2 = flint_malloc(sizeof(slong) * m);

        /* Generate matrix with known bounds on rank */
        if (n_randint(state, 2))
        {
            GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));

            if (lu_impl != ((gr_method_mat_lu_op) gr_mat_lu_classical) &&
                GR_SUCCESS == gr_mat_lu_classical(&rank_lower_bound, P2, LU2, A, 0, ctx))
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
            _gr_mat_randrank(A, state, r, 5, ctx);
            if (n_randint(state, 2))
            {
                d = n_randint(state, 2*m*n + 1);
                GR_MUST_SUCCEED(gr_mat_randops(A, state, d, ctx));
            }

            if (gr_ctx_is_finite_characteristic(ctx) == T_FALSE)
                rank_lower_bound = rank_upper_bound = r;
            else
            {
                rank_lower_bound = 0;
                rank_upper_bound = r;
            }
        }

        status = lu_impl(&rank, P, LU, A, 0, ctx);

        if (status == GR_SUCCESS)
        {
            /* Check shape of solution */
            check(P, LU, A, rank, ctx);

            /* Check rank */
            if (rank < rank_lower_bound || rank > rank_upper_bound)
            {
                flint_printf("FAIL\n");
                gr_ctx_println(ctx);
                flint_printf("wrong rank!\n");
                flint_printf("rank = %wd\n", rank);
                flint_printf("rank bounds = [%wd, %wd]\n", rank_lower_bound, rank_upper_bound);
                flint_printf("\n\nA:");
                gr_mat_print(A, ctx);
                flint_printf("\n\nLU:");
                gr_mat_print(LU, ctx);
                flint_abort();
            }
        }

        /* ----------------------------------------------------- */
        /* Rank check for square matrices                        */
        /* ----------------------------------------------------- */
        if (m == n)
        {
            status = lu_impl(&rank, P, LU, A, 1, ctx);

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
                    gr_ctx_println(ctx);
                    flint_printf("rank check\n");
                    flint_printf("rank = %wd\n", rank);
                    flint_printf("rank bounds = [%wd, %wd]\n", rank_lower_bound, rank_upper_bound);
                    flint_printf("\n\nA:\n");
                    gr_mat_print(A, ctx);
                    flint_printf("\n\nLU:\n");
                    gr_mat_print(LU, ctx);
                    flint_abort();
                }

                if (rank != 0)
                    check(P, LU, A, rank, ctx);
            }
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(LU, ctx);
        gr_mat_clear(LU2, ctx);
        flint_free(P);
        flint_free(P2);

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}
