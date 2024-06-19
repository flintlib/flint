/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

void _gr_poly_test_inv_series(gr_method_poly_unary_trunc_op inv_series_impl,
    flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ptr A, invA, AinvA, one;
        slong lenA, len;
        int status = GR_SUCCESS;
        gr_ctx_t my_ctx;
        gr_ctx_struct * ctx;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        lenA = 1 + n_randint(state, maxn);
        len = 1 + n_randint(state, maxn);

        A = gr_heap_init_vec(lenA, ctx);
        invA = gr_heap_init_vec(len, ctx);
        AinvA = gr_heap_init_vec(len, ctx);
        one = gr_heap_init_vec(len, ctx);

        GR_MUST_SUCCEED(_gr_vec_randtest(A, state, lenA, ctx));
        GR_MUST_SUCCEED(_gr_vec_randtest(invA, state, len, ctx));

        status |= inv_series_impl(invA, A, lenA, len, ctx);

        if (status == GR_SUCCESS)
        {
            status |= _gr_poly_mullow(AinvA, A, lenA, invA, len, len, ctx);

            status |= _gr_vec_zero(one, len, ctx);
            status |= gr_one(one, ctx);

            if (status == GR_SUCCESS && _gr_vec_equal(AinvA, one, len, ctx) == T_FALSE)
            {
                flint_printf("FAIL:\n");
                gr_ctx_println(ctx);
                flint_printf("lenA = %wd, len = %wd\n\n", lenA, len);
                flint_printf("A:\n"); _gr_vec_print(A, lenA, ctx); flint_printf("\n\n");
                flint_printf("invA:\n"); _gr_vec_print(invA, len, ctx); flint_printf("\n\n");
                flint_printf("AinvA:\n"); _gr_vec_print(AinvA, len, ctx); flint_printf("\n\n");
                flint_abort();
            }
        }

        gr_heap_clear_vec(A, lenA, ctx);
        gr_heap_clear_vec(invA, len, ctx);
        gr_heap_clear_vec(AinvA, len, ctx);
        gr_heap_clear_vec(one, len, ctx);

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}
