/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"
#include "ca_mat.h"

/* Defined in t-rref.c, t-rref_fflu.c and t-rref_lu.c */
#ifndef ca_mat_randrowops
#define ca_mat_randrowops ca_mat_randrowops
void
ca_mat_randrowops(ca_mat_t mat, flint_rand_t state, slong count, ca_ctx_t ctx)
{
    slong c, i, j, k;
    slong m = mat->r;
    slong n = mat->c;

    if (mat->r == 0 || mat->c == 0)
        return;

    for (c = 0; c < count; c++)
    {
        if ((i = n_randint(state, m)) == (j = n_randint(state, m)))
            continue;
        if (n_randint(state, 2))
            for (k = 0; k < n; k++)
                ca_add(ca_mat_entry(mat, j, k), ca_mat_entry(mat, j, k), ca_mat_entry(mat, i, k), ctx);
        else
            for (k = 0; k < n; k++)
                ca_sub(ca_mat_entry(mat, j, k), ca_mat_entry(mat, j, k), ca_mat_entry(mat, i, k), ctx);
    }
}
#endif

/* On success, R must be in reduced row echelon form: in each of the first
   rank rows, the first entry that is not provably zero is provably one and
   is the only entry not provably zero in its column; the remaining rows are
   zero. This catches returning success on an unfinished computation. */
static int
_ca_mat_check_rref_shape(const ca_mat_t R, slong rank, ca_ctx_t ctx)
{
    slong i, j, k, prev;

    prev = -1;

    for (i = 0; i < ca_mat_nrows(R); i++)
    {
        for (j = 0; j < ca_mat_ncols(R); j++)
            if (ca_check_is_zero(ca_mat_entry(R, i, j), ctx) != T_TRUE)
                break;

        if (i >= rank)
        {
            if (j != ca_mat_ncols(R))
                return 0;
            continue;
        }

        if (j == ca_mat_ncols(R) || j <= prev)
            return 0;

        if (ca_check_is_one(ca_mat_entry(R, i, j), ctx) != T_TRUE)
            return 0;

        for (k = 0; k < ca_mat_nrows(R); k++)
            if (k != i && ca_check_is_zero(ca_mat_entry(R, k, j), ctx) != T_TRUE)
                return 0;

        prev = j;
    }

    return 1;
}

TEST_FUNCTION_START(ca_mat_rref_lu, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, R, R2;
        fmpq_mat_t AQ, RQ;
        slong rank, rankq, r, c;
        int success;

        ca_ctx_init(ctx);

        r = n_randint(state, 10);
        c = n_randint(state, 10);

        ca_mat_init(A, r, c, ctx);
        ca_mat_init(R, r, c, ctx);
        ca_mat_init(R2, r, c, ctx);
        fmpq_mat_init(AQ, r, c);
        fmpq_mat_init(RQ, r, c);

        fmpq_mat_randtest(AQ, state, 10);

        ca_mat_set_fmpq_mat(A, AQ, ctx);

        rankq = fmpq_mat_rref(RQ, AQ);

        if (n_randint(state, 2))
        {
            ca_mat_set(R, A, ctx);
            success = ca_mat_rref_lu(&rank, R, R, ctx);
        }
        else
        {
            success = ca_mat_rref_lu(&rank, R, A, ctx);
        }

        if (success)
        {
            ca_mat_set_fmpq_mat(R2, RQ, ctx);

            if (rank != rankq || ca_mat_check_equal(R, R2, ctx) == T_FALSE)
            {
                flint_printf("FAIL:\n");
                flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                flint_printf("R: "); ca_mat_print(R, ctx); flint_printf("\n");
                flint_printf("R2: "); ca_mat_print(R2, ctx); flint_printf("\n");
                flint_printf("rank = %wd, %wd\n\n", rank, rankq);
                flint_abort();
            }
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(R, ctx);
        ca_mat_clear(R2, ctx);
        fmpq_mat_clear(AQ);
        fmpq_mat_clear(RQ);

        ca_ctx_clear(ctx);
    }

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_mat_t A, B, R, R2;
        slong rank1, rank2, r, c;
        int success;

        ca_ctx_init(ctx);

        r = n_randint(state, 5);
        c = n_randint(state, 5);

        ca_mat_init(A, r, c, ctx);
        ca_mat_init(B, r, c, ctx);
        ca_mat_init(R, r, c, ctx);
        ca_mat_init(R2, r, c, ctx);

        ca_mat_randtest(A, state, 1, 5, ctx);
        ca_mat_set(B, A, ctx);
        ca_mat_randrowops(B, state, 1 + n_randint(state, 20), ctx);

        success = ca_mat_rref_lu(&rank1, R, A, ctx);

        if (success && !_ca_mat_check_rref_shape(R, rank1, ctx))
        {
            flint_printf("FAIL (not in rref form):\n");
            flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
            flint_printf("R: "); ca_mat_print(R, ctx); flint_printf("\n");
            flint_printf("rank = %wd\n\n", rank1);
            flint_abort();
        }

        if (success)
        {
            success = ca_mat_rref_lu(&rank2, R2, B, ctx);

            if (success)
            {
                if (rank1 != rank2 || ca_mat_check_equal(R, R2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("A: "); ca_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("B: "); ca_mat_print(B, ctx); flint_printf("\n");
                    flint_printf("R: "); ca_mat_print(R, ctx); flint_printf("\n");
                    flint_printf("R2: "); ca_mat_print(R2, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank1, rank2);
                    flint_abort();
                }
            }
        }

        ca_mat_clear(A, ctx);
        ca_mat_clear(B, ctx);
        ca_mat_clear(R, ctx);
        ca_mat_clear(R2, ctx);

        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
