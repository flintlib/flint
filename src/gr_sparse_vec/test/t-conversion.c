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

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }

int
test_conversion(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    slong N = 100;
    slong n_tests = 20;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_sparse_vec_t vec, vec2;
    gr_vec_t dvec, dvec2;

    gr_vec_init(dvec, N, ctx);
    gr_vec_init(dvec2, N, ctx);
    gr_sparse_vec_init(vec, N, ctx);
    gr_sparse_vec_init(vec2, N, ctx);

    for (i = 0; i < n_tests; i++)
    {
        status |= _gr_vec_randtest(GR_VEC_ENTRY(dvec, 0, sz), state, N, ctx);
        status |= gr_sparse_vec_set_vec(vec, GR_VEC_ENTRY(dvec, 0, sz), N, ctx);
        status |= gr_vec_set_sparse_vec(GR_VEC_ENTRY(dvec2, 0, sz), vec, ctx);
        if (T_TRUE != _gr_vec_equal(GR_VEC_ENTRY(dvec, 0, sz), GR_VEC_ENTRY(dvec2, 0, sz), N, ctx))
            return GR_TEST_FAIL;
    }
    gr_sparse_vec_clear(vec, ctx);
    gr_vec_clear(dvec, ctx);
    gr_vec_clear(dvec2, ctx);
    return status;
}

TEST_FUNCTION_START(gr_sparse_vec_conversion, state)
{   
    gr_ctx_t ctx;
    gr_ctx_init_random(ctx, state);
    CHECK_TEST(test_conversion(state, ctx), "Conversion from and to dense");
    gr_ctx_clear(ctx);
    TEST_FUNCTION_END(state);
}
