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

void _gr_poly_test_div_series(gr_method_poly_binary_trunc_op div_series_impl,
    flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ptr A, B, AdivB, A2;
        slong lenA, lenB, len;
        int status = GR_SUCCESS;
        gr_ctx_t ctx2;
        gr_ctx_struct * ctxptr;

        if (ctx == NULL)
        {
            gr_ctx_init_random(ctx2, state);
            ctxptr = ctx2;
        }
        else
            ctxptr = ctx;

        lenA = 1 + n_randint(state, maxn);
        lenB = 1 + n_randint(state, maxn);
        len = 1 + n_randint(state, maxn);

        A = gr_heap_init_vec(lenA, ctxptr);
        B = gr_heap_init_vec(lenB, ctxptr);
        AdivB = gr_heap_init_vec(len, ctxptr);
        A2 = gr_heap_init_vec(len, ctxptr);

        GR_MUST_SUCCEED(_gr_vec_randtest(A, state, lenA, ctxptr));
        GR_MUST_SUCCEED(_gr_vec_randtest(B, state, lenB, ctxptr));
        GR_MUST_SUCCEED(_gr_vec_randtest(AdivB, state, len, ctxptr));

        status |= div_series_impl(AdivB, A, lenA, B, lenB, len, ctxptr);

        if (status == GR_SUCCESS)
        {
            status |= _gr_poly_mullow(A2, AdivB, len, B, lenB, len, ctxptr);

            if (status == GR_SUCCESS && _gr_poly_equal(A2, len, A, FLINT_MIN(lenA, len), ctxptr) == T_FALSE)
            {
                flint_printf("FAIL:\n");
                gr_ctx_println(ctxptr);
                flint_printf("lenA = %wd, len = %wd\n\n", lenA, len);
                flint_printf("A:\n"); _gr_vec_print(A, lenA, ctxptr); flint_printf("\n\n");
                flint_printf("B:\n"); _gr_vec_print(B, lenB, ctxptr); flint_printf("\n\n");
                flint_printf("AdivB:\n"); _gr_vec_print(AdivB, len, ctxptr); flint_printf("\n\n");
                flint_printf("A2:\n"); _gr_vec_print(A2, len, ctxptr); flint_printf("\n\n");
                flint_abort();
            }
        }

        gr_heap_clear_vec(A, lenA, ctxptr);
        gr_heap_clear_vec(B, len, ctxptr);
        gr_heap_clear_vec(AdivB, len, ctxptr);
        gr_heap_clear_vec(A2, len, ctxptr);

        if (ctx == NULL)
            gr_ctx_clear(ctxptr);
    }
}
