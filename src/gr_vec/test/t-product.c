/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_vec.h"

int
test_product(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    gr_ptr vec, x1, x2;
    slong len;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    if (gr_ctx_is_finite(ctx) == T_TRUE || gr_ctx_has_real_prec(ctx) == T_TRUE)
        len = n_randint(state, 120);
    else
        len = n_randint(state, 20);

    GR_TMP_INIT_VEC(vec, len, ctx);
    GR_TMP_INIT2(x1, x2, ctx);

    status |= _gr_vec_randtest(vec, state, len, ctx);
    status |= _gr_vec_product_serial(x1, vec, len, ctx);

    if (gr_ctx_is_threadsafe(ctx) == T_TRUE && n_randint(state, 4) == 0)
        flint_set_num_threads(1 + n_randint(state, 4));
    else
        flint_set_num_threads(1);

    switch (which)
    {
        case 0: status |= _gr_vec_product_bsplit(x2, vec, len, 2 + n_randint(state, 8), ctx); break;
        case 1: status |= _gr_vec_product_bsplit_parallel(x2, vec, len, 2 + n_randint(state, 8), ctx); break;
        case 2: status |= _gr_vec_product_parallel(x2, vec, len, ctx); break;
        default: status |= _gr_vec_product(x2, vec, len, ctx); break;
    }

    if (status == GR_SUCCESS && gr_equal(x1, x2, ctx) == T_FALSE)
    {
        flint_printf("FAIL\n");
        printf("which = %d\n", which);
        printf("vec = "); _gr_vec_print(vec, len, ctx); printf("\n");
        printf("x1 = "); gr_println(x1, ctx); printf("\n");
        printf("x2 = "); gr_println(x2, ctx); printf("\n");
        flint_abort();
    }

    GR_TMP_CLEAR2(x1, x2, ctx);
    GR_TMP_CLEAR_VEC(vec, len, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_vec_product, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test_product(state, 0);
        test_product(state, 1);
        test_product(state, 2);
        test_product(state, 3);
    }

    TEST_FUNCTION_END(state);
}
