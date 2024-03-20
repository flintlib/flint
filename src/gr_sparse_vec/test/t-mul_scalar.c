/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_sparse_vec.h"
#include "fmpz.h"
#include "fmpq.h"

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }


int
test_mul_scalar(flint_rand_t state, gr_ctx_t ctx)
{
    slong i, c;
    ulong uc;
    slong N = 100;
    slong n_tests = 20;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    gr_sparse_vec_t vec, vec2;
    gr_ptr dvec, dvec2;
    gr_ptr temp;
    fmpz_t zc;
    fmpq_t qc;

    GR_TMP_INIT(temp, ctx);
    fmpz_init(zc);
    fmpq_init(qc);
    dvec = flint_malloc(N * sz);
    dvec2 = flint_malloc(N * sz);
    _gr_vec_init(dvec, N, ctx);
    _gr_vec_init(dvec2, N, ctx);
    gr_sparse_vec_init(vec, N, ctx);
    gr_sparse_vec_init(vec2, N, ctx);

    for (i = 0; i < 6 * n_tests; i++)
    {
        if (i % 6 == 4 && gr_ctx_is_field(ctx) != T_TRUE)
            continue;
        status |= gr_sparse_vec_randtest(vec, 10, 0, state, ctx);
        status |= gr_vec_set_sparse_vec(dvec, vec, ctx);

        switch(i % 6)
        {
        case 0: 
            //flint_printf("Testing scalar\n");
            status |= gr_randtest_not_zero(temp, state, ctx);
            status |= gr_sparse_vec_mul_scalar(vec2, vec, temp, ctx);
            status |= _gr_vec_mul_scalar(dvec2, dvec, N, temp, ctx);
            break;
        case 1:
            //flint_printf("Testing scalar_si\n");
            c = n_randint(state, 0);
            status |= gr_sparse_vec_mul_scalar_si(vec2, vec, c, ctx);
            status |= _gr_vec_mul_scalar_si(dvec2, dvec, N, c, ctx);
            break;
        case 2:
            //flint_printf("Testing scalar_ui\n");
            uc = n_randint(state, 0);
            status |= gr_sparse_vec_mul_scalar_ui(vec2, vec, uc, ctx);
            status |= _gr_vec_mul_scalar_ui(dvec2, dvec, N, uc, ctx);
            break;
        case 3:
            //flint_printf("Testing scalar_fmpz\n");
            fmpz_randtest_not_zero(zc, state, 32);
            status |= gr_sparse_vec_mul_scalar_fmpz(vec2, vec, zc, ctx);
            status |= _gr_vec_mul_scalar_fmpz(dvec2, dvec, N, zc, ctx);
            break;
        case 4:
            //flint_printf("Testing scalar_fmpq\n");
            fmpq_randtest_not_zero(qc, state, 32);
            //fmpq_print(qc); flint_printf("\n");
            status |= gr_sparse_vec_mul_scalar_fmpq(vec2, vec, qc, ctx);
            //_gr_vec_print(dvec, N, ctx); flint_printf("\n");
            status |= _gr_vec_mul_scalar_fmpq(dvec2, dvec, N, qc, ctx);
            //_gr_vec_print(dvec2, N, ctx); flint_printf("\n\n\n");
           break;
        case 5:
            //flint_printf("Testing scalar_2exp_si\n");
            c = n_randint(state, 32) + 1;
            status |= gr_sparse_vec_mul_scalar_2exp_si(vec2, vec, c, ctx);
            status |= _gr_vec_mul_scalar_2exp_si(dvec2, dvec, N, c, ctx);
            break;
        }
        //gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
        //gr_sparse_vec_print_nz(vec2, ctx); flint_printf("\n");
 
        status |= gr_vec_set_sparse_vec(dvec, vec2, ctx);
        if (T_FALSE == _gr_vec_equal(dvec, dvec2, N, ctx))
        {
            _gr_vec_print(dvec, N, ctx); flint_printf("\n");
            _gr_vec_print(dvec2, N, ctx); flint_printf("\n");
            status |= gr_sparse_vec_set_vec(vec, dvec2, N, ctx);
            gr_sparse_vec_print_nz(vec, ctx); flint_printf("\n");
            gr_sparse_vec_print_nz(vec2, ctx); flint_printf("\n");
            return GR_TEST_FAIL;
        }
    }

    gr_sparse_vec_clear(vec, ctx);
    gr_sparse_vec_clear(vec2, ctx);
    _gr_vec_clear(dvec, N, ctx);
    _gr_vec_clear(dvec2, N, ctx);
    flint_free(dvec);
    flint_free(dvec2);
    GR_TMP_CLEAR(temp, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_vec_mul_div_scalar, state)
{   
    slong i;
    gr_ctx_t ctx;
    for (i = 0; i < 16; ++i)
    {
        while (1)
        {
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_zero_ring(ctx) == T_TRUE)
                gr_ctx_clear(ctx);
            else
                break;
        }
        //gr_ctx_println(ctx);
        CHECK_TEST(test_mul_scalar(state, ctx), "Scalar multiplication and division");
        gr_ctx_clear(ctx);
    }
    TEST_FUNCTION_END(state);
}
