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

int test_solve(flint_rand_t state, gr_ctx_t ctx)
{
    int iter, ret;
    slong rep, nreps = 100, r, c, i;
    gr_coo_mat_t Aorig;
    gr_lil_mat_t A, At;
    gr_ptr x, x2, b, Atb, Ax, AtAx;
    gr_mat_t dmat;
    slong niters[6] = {0, 0, 0, 0, 0, 0};
    slong psol[6] = {0, 0, 0, 0, 0, 0};
    slong nosol[6] = {0, 0, 0, 0, 0, 0};
    slong sol[6] = {0, 0, 0, 0, 0, 0};
    /* double elapsed[6] = {0, 0, 0, 0, 0}; */
    char *names[6] = {"rref", "lu", "Lanczos", "block Lanczos", "Wiedemann", "block Wiedemann"};
    int status = GR_SUCCESS;
    /* struct timeval start, end; */
    
    flint_printf("solving Ax = b....");
    fflush(stdout);
    
    for (rep = 0; rep < nreps; rep++)
    {
        if (rep % 5==0) {flint_printf("."); fflush(stdout);}

        c = r = 50 + n_randint(state, 5);

        gr_coo_mat_init(Aorig, r, c, ctx);
        gr_lil_mat_init(A, r, c, ctx);
        gr_mat_init(dmat, r, c, ctx);
        GR_TMP_INIT_VEC(x, c, ctx);
        GR_TMP_INIT_VEC(x2, c, ctx);
        GR_TMP_INIT_VEC(b, r, ctx);
        GR_TMP_INIT_VEC(Ax, r, ctx);

        // Set up a solvable problem
        status |= gr_coo_mat_randtest_prob(Aorig, .2, state, ctx);
        status |= gr_lil_mat_set_coo_mat(A, Aorig, ctx);
        status |= _gr_vec_randtest(x, state, c, ctx);
        status |= gr_lil_mat_mul_vec(b, A, x, ctx);
        status |= gr_mat_set_lil_mat(dmat, A, ctx);
        // flint_printf("A = "); gr_mat_print(dmat, ctx); flint_printf("\n");
        // flint_printf("x = "); _gr_vec_print(x, c, ctx); flint_printf("\n");
        // flint_printf("b = "); _gr_vec_print(b, r, ctx); flint_printf("\n");

        for (i = 0; i < 6; ++i)
        {
            iter = 0;
            /* gettimeofday(&start, NULL); */
            // TODO: rref and lu solving
            if (i == 0 || i == 1)
                continue;
            switch (i) 
            {
            //case 0: ret = gr_lil_mat_solve_rref(x2, A, b, ctx); break;
            //case 1: ret = gr_lil_mat_solve_lu(x2, A, b, ctx); break;
            case 2: do ret = gr_lil_mat_solve_lanczos(x2, A, b, state, ctx); while (ret == GR_UNABLE && ++iter < 3); break;
            case 3: do ret = gr_lil_mat_solve_block_lanczos(x2, A, b, 8, state, ctx); while (ret == GR_UNABLE && ++iter < 3); break;
            case 4: ret = gr_lil_mat_solve_wiedemann(x2, A, b, ctx); break;
            case 5: do ret = gr_lil_mat_solve_block_wiedemann(x2, A, b, 8, state, ctx); while (ret == GR_TEST_FAIL && ++iter < 3); break;
            }
            // /* gettimeofday(&end, NULL);
            // elapsed[i] += (end.tv_sec - start.tv_sec) + .000001*(end.tv_usec-start.tv_usec); */
            if (ret == GR_UNABLE) nosol[i] += 1;
            else
            {
                niters[i] += iter;
                status |= gr_lil_mat_mul_vec(Ax, A, x2, ctx);
                if (_gr_vec_equal(b, Ax, A->r, ctx) == T_FALSE) 
                {
                    if (i == 2 || i == 3)
                    {
                        gr_lil_mat_init(At, c, r, ctx);
                        GR_TMP_INIT_VEC(AtAx, c, ctx);
                        GR_TMP_INIT_VEC(Atb, c, ctx);
                        status |= gr_lil_mat_transpose(At, A, ctx);
                        status |= gr_lil_mat_mul_vec(AtAx, At, Ax, ctx);
                        status |= gr_lil_mat_mul_vec(Atb, At, b, ctx);
                        if (_gr_vec_equal(AtAx, Atb, A->c, ctx) != T_TRUE)
                        {
                            flint_printf("FAIL on %s: AtAx != Atb\n", names[i]);
                            abort();
                        } 
                        else psol[i] += 1;
                        GR_TMP_CLEAR_VEC(AtAx, c, ctx);
                        GR_TMP_CLEAR_VEC(Atb, c, ctx);
                        gr_lil_mat_clear(At, ctx);
                    }
                    else
                    {
                        flint_printf("FAIL on %s: Ax != b\n", names[i]);
                        return GR_TEST_FAIL;
                    }
                }
                else
                {
                    sol[i] += 1;
                }
            }
        }

        GR_TMP_CLEAR_VEC(x, c, ctx);
        GR_TMP_CLEAR_VEC(x2, c, ctx);
        GR_TMP_CLEAR_VEC(b, r, ctx);
        GR_TMP_CLEAR_VEC(Ax, r, ctx);
        gr_lil_mat_clear(A, ctx);
        gr_mat_clear(dmat, ctx);
    }
    
    flint_printf("PASS\n");
    for (i = 0; i < 6; ++i)
    {
        flint_printf("Solved %d with %s\n", sol[i], names[i]);
        /* flint_printf("\tAverage time: %lf\n", elapsed[i]/nreps); */
        if (nosol[i])
            flint_printf("\tFound no solution for %wd/%wd examples\n", nosol[i], nreps);
        if (psol[i])    
            flint_printf("\tFound pseudo-solution for %wd/%wd examples\n", psol[i], nreps);
        if (niters[i])
            flint_printf("\tRequired %f extra iters per solution (on average).\n", (double)niters[i]/nreps);
    }
    return status;
}


TEST_FUNCTION_START(gr_sparse_mat_solve, state)
{   
    slong i;
    gr_ctx_t ctx;
    for (i = 0; i < 1; ++i)
    {
        while (1)
        {
            //gr_ctx_init_fmpz(ctx);
            gr_ctx_init_random(ctx, state);
            gr_ctx_init_fq_nmod(ctx, 65521, 1, "a");
            if (gr_ctx_is_zero_ring(ctx) == T_TRUE || gr_ctx_is_exact(ctx) != T_TRUE || gr_ctx_is_field(ctx) != T_TRUE || gr_ctx_is_finite(ctx) != T_TRUE)
                gr_ctx_clear(ctx);
            else
                break;
        }
        CHECK_TEST(test_solve(state, ctx), "Sparse matrix solving");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}

