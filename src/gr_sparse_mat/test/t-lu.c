/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "test_helpers.h"
#include "gr_sparse_mat.h"
/* #include <sys/time.h> */

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }

int test_lu(flint_rand_t state, gr_ctx_t ctx)
{
    slong rep, r, c, i, j, rk, *P, *Q;
    gr_ptr val;
    gr_coo_mat_t Aorig;
    gr_lil_mat_t A, LU, L, U;
    gr_mat_t dL, dU, dLU;
    int status = GR_SUCCESS;
    
    flint_printf("decomposing PAQ = LU....");
    fflush(stdout);

    for (rep = 0; rep < 200; rep++)
    {
        if (rep % 20 == 0) {flint_printf("."); fflush(stdout);}

        r = n_randint(state, 20);
        c = n_randint(state, 20);

        P = flint_malloc(r*sizeof(*P));
        Q = flint_malloc(c*sizeof(*P));

        gr_coo_mat_init(Aorig, r, c, ctx);
        gr_lil_mat_init(A, r, c, ctx);
        gr_lil_mat_init(LU, r, c, ctx);
        gr_lil_mat_init(L, r, c, ctx);
        gr_lil_mat_init(U, r, c, ctx);

        status |= gr_coo_mat_randtest(Aorig, FLINT_MIN(r * 3, c * 3), 0, T_TRUE, state, ctx);
        status |= gr_lil_mat_set_coo_mat(A, Aorig, ctx);
        status |= gr_lil_mat_lu(&rk, P, Q, L, U, A, ctx);
        
        /* Check that L is lower triangular (with ones on diagonal up to rank) */
        for (i = 0; i < r; ++i) 
        {
            val = gr_sparse_vec_find_entry(&L->rows[i], i, ctx);
            if (i < rk && (val == NULL || gr_is_one(val, ctx) == T_FALSE))
            {
                flint_printf("FAIL: L does not have unit diagonal up to the rank\n");
                 status |= gr_lil_mat_print_nz(L, ctx); flint_printf("\n");
                 return GR_TEST_FAIL;
            }
            for (j = 0; j < L->rows[i].nnz; ++j) 
            {
                if (L->rows[i].inds[j] > i) 
                {
                    flint_printf("FAIL: L not lower triangular\n");
                    status |= gr_lil_mat_print_nz(L, ctx); flint_printf("\n");
                    return GR_TEST_FAIL;
                }
                if (L->rows[i].inds[j] >= rk) 
                {
                    flint_printf("FAIL: L not trivial past the rank\n");
                    status |= gr_lil_mat_print_nz(L, ctx); flint_printf("\n");
                    /*gr_lil_mat_print_pretty(L, ctx);*/
                    return GR_TEST_FAIL;
                }
            }
        }
        /* Check that U is upper triangular (with nonzero diagonal up to rank) */
        for (i = 0; i < r; ++i) 
        {
            val = gr_sparse_vec_find_entry(&U->rows[i], i, ctx);
            if (i < rk && (val == NULL || gr_is_zero(val, ctx) == T_TRUE))
            {
                flint_printf("FAIL: U does not have nonzero diagonal\n");
                status |= gr_lil_mat_print_nz(U, ctx); flint_printf("\n");
                return GR_TEST_FAIL;
            }
            if (i >= rk && U->rows[i].nnz != 0) 
            {
                flint_printf("FAIL: U not trivial past the rank\n");
                status |= gr_lil_mat_print_nz(U, ctx); flint_printf("\n");
                return GR_TEST_FAIL;
            }
            for (j = 0; j < U->rows[i].nnz; ++j) 
            {
                if (U->rows[i].inds[j] < i) 
                {
                    flint_printf("FAIL: U not upper triangular\n");
                    status |= gr_lil_mat_print_nz(U, ctx); flint_printf("\n");
                    return GR_TEST_FAIL;
                }
            }
        }

        //flint_printf("L = "); status |= gr_lil_mat_print_nz(L, ctx); flint_printf("\n");
        //flint_printf("U = "); status |= gr_lil_mat_print_nz(U, ctx); flint_printf("\n");
        
        // Check that PAQ = LU
        gr_mat_init(dL, r, c, ctx);
        gr_mat_init(dU, r, c, ctx);
        gr_mat_init(dLU, r, c, ctx);
        status |= gr_mat_set_lil_mat(dL, L, ctx);
        status |= gr_mat_set_lil_mat(dU, U, ctx);
        status |= gr_mat_mul(dLU, dL, dU, ctx);
        status |= gr_lil_mat_set_mat(LU, dLU, ctx);
        status |= gr_lil_mat_permute_rows(A, P, ctx);
        status |= gr_lil_mat_permute_cols(A, Q, ctx);
        if (status != GR_SUCCESS || gr_lil_mat_equal(A, LU, ctx) == T_FALSE) 
        {
            flint_printf("FAIL: PAQ != LU\n");
            flint_printf("PAQ = "); status |= gr_lil_mat_print_nz(A, ctx); flint_printf("\n");
            flint_printf("LU = "); status |= gr_lil_mat_print_nz(LU, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }

        flint_free(P);
        flint_free(Q);
        gr_coo_mat_clear(Aorig, ctx);
        gr_lil_mat_clear(A, ctx);
        gr_lil_mat_clear(U, ctx);
        gr_lil_mat_clear(L, ctx);
        gr_lil_mat_clear(LU, ctx);
        gr_mat_clear(dL, ctx);
        gr_mat_clear(dU, ctx);
        gr_mat_clear(dLU, ctx);
    }
    
    flint_printf("PASS\n");
    return status;
}


TEST_FUNCTION_START(gr_sparse_mat_lu, state)
{   
    slong i;
    gr_ctx_t ctx;
    for (i = 0; i < 1; ++i)
    {
        while (1)
        {
            //gr_ctx_init_fmpz(ctx);
            //gr_ctx_init_random(ctx, state);
            gr_ctx_init_fq_nmod(ctx, 65521, 1, "a");
            if (gr_ctx_is_zero_ring(ctx) == T_TRUE || gr_ctx_is_exact(ctx) != T_TRUE || gr_ctx_is_field(ctx) != T_TRUE || gr_ctx_is_finite(ctx) != T_TRUE)
                gr_ctx_clear(ctx);
            else
                break;
        }
        CHECK_TEST(test_lu(state, ctx), "Sparse matrix LU decomposition");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
