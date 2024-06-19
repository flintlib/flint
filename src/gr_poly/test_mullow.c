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

void _gr_poly_test_mullow(gr_method_poly_binary_trunc_op mullow_impl, gr_method_poly_binary_trunc_op mullow_ref,
    flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    if (mullow_ref == NULL)
        mullow_ref = (gr_method_poly_binary_trunc_op) _gr_poly_mullow_generic;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ptr A, B, C, D;
        slong n1, n2, n;
        int status = GR_SUCCESS;
        gr_ctx_t my_ctx;
        gr_ctx_struct * ctx;
        int squaring;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        squaring = n_randint(state, 2);

        n1 = 1 + n_randint(state, 1 + maxn);
        n2 = squaring ? n1 : 1 + n_randint(state, 1 + maxn);
        n = 1 + n_randint(state, 1 + maxn);
        n = FLINT_MIN(n, n1 + n2 - 1);

        A = gr_heap_init_vec(n1, ctx);
        B = squaring ? A : gr_heap_init_vec(n2, ctx);
        C = gr_heap_init_vec(n, ctx);
        D = gr_heap_init_vec(n, ctx);

        GR_MUST_SUCCEED(_gr_vec_randtest(A, state, n1, ctx));
        if (!squaring)
            GR_MUST_SUCCEED(_gr_vec_randtest(B, state, n2, ctx));
        GR_MUST_SUCCEED(_gr_vec_randtest(C, state, n, ctx));

#if 0
        /* Useful for debugging */
        slong i;
        for (i = 0; i < n1; i++)
            GR_IGNORE(gr_set_ui(GR_ENTRY(A, i, ctx->sizeof_elem), i + 1, ctx));
        for (i = 0; i < n2; i++)
            GR_IGNORE(gr_set_ui(GR_ENTRY(B, i, ctx->sizeof_elem), i + 1, ctx));
#endif

        status |= mullow_impl(C, A, n1, B, n2, n, ctx);

        if (status == GR_SUCCESS)
            status |= mullow_ref(D, A, n1, B, n2, n, ctx);

        if (status == GR_SUCCESS && _gr_vec_equal(C, D, n, ctx) == T_FALSE)
        {
            flint_printf("FAIL:\n");
            gr_ctx_println(ctx);
            flint_printf("len1 = %wd, len2 = %wd, len = %wd, squaring = %d\n\n", n1, n2, n, squaring);
            flint_printf("A:\n"); _gr_vec_print(A, n1, ctx); flint_printf("\n\n");
            flint_printf("B:\n"); _gr_vec_print(B, n2, ctx); flint_printf("\n\n");
            flint_printf("C:\n"); _gr_vec_print(C, n, ctx); flint_printf("\n\n");
            flint_printf("D:\n"); _gr_vec_print(D, n, ctx); flint_printf("\n\n");
            flint_abort();
        }

        gr_heap_clear_vec(A, n1, ctx);
        if (!squaring)
            gr_heap_clear_vec(B, n2, ctx);
        gr_heap_clear_vec(C, n, ctx);
        gr_heap_clear_vec(D, n, ctx);

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}
