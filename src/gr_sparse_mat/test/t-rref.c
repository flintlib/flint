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

int test_rref(flint_rand_t state, gr_ctx_t ctx)
{
    slong rep, r, c, sparse_rank, dense_rank;
    gr_coo_mat_t Aorig;
    gr_lil_mat_t A, B, R;
    gr_mat_t dA, dR;
    int status = GR_SUCCESS;

    flint_printf("converting A to reduced row echelon form....");
    fflush(stdout);
    
    for (rep = 0; rep < 200; rep++)
    {
        if (rep % 20 == 0) {flint_printf("."); fflush(stdout);}

        r = n_randint(state, 20);
        c = n_randint(state, 20);
        gr_coo_mat_init(Aorig, r, c, ctx);
        gr_lil_mat_init(A, r, c, ctx);
        gr_lil_mat_init(B, r, c, ctx);
        gr_lil_mat_init(R, r, c, ctx);
        gr_mat_init(dA, r, c, ctx);
        gr_mat_init(dR, r, c, ctx);

        status |= gr_coo_mat_randtest(Aorig, FLINT_MIN(r * 3, c * 3), 0, T_TRUE, state, ctx);
        if (status != GR_SUCCESS) 
        {
            flint_printf("Some failure!\n");
            return GR_TEST_FAIL;
        }
        status |= gr_lil_mat_set_coo_mat(A, Aorig, ctx);
        status |= gr_mat_set_lil_mat(dA, A, ctx);

        status |= gr_lil_mat_rref(&sparse_rank, R, A, ctx);
        status |= gr_mat_rref(&dense_rank, dR, dA, ctx);
        status |= gr_lil_mat_set_mat(B, dR, ctx);

        if (status != GR_SUCCESS || gr_lil_mat_equal(R, B, ctx) == T_FALSE) 
        {
            flint_printf("FAIL!\n");
            flint_printf("A = "); status |= gr_lil_mat_print_nz(A, ctx); flint_printf("\n");
            flint_printf("R = "); status |= gr_lil_mat_print_nz(R, ctx); flint_printf("\n");
            flint_printf("dR = "); status |= gr_lil_mat_print_nz(B, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }

        gr_coo_mat_clear(Aorig, ctx);
        gr_lil_mat_clear(A, ctx);
        gr_lil_mat_clear(B, ctx);
        gr_lil_mat_clear(R, ctx);
        gr_mat_clear(dA, ctx);
        gr_mat_clear(dR, ctx);
    }
    
    flint_printf("PASS\n");
    return status;
}


TEST_FUNCTION_START(gr_sparse_mat_rref, state)
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
        CHECK_TEST(test_rref(state, ctx), "Sparse matrix reduced row echelon form");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}

