/*
    Copyright (C) 2022, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you grn redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

void gr_mat_test_charpoly(gr_method_mat_unary_op_get_scalar charpoly_impl, flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t my_ctx;
        gr_ctx_ptr ctx;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        {
            gr_mat_t A;
            gr_ptr cp, cp2;
            slong n;
            int status = GR_SUCCESS;

            n = n_randint(state, maxn + 1);

            gr_mat_init(A, n, n, ctx);
            cp = gr_heap_init_vec(n + 1, ctx);
            cp2 = gr_heap_init_vec(n + 1, ctx);

            status |= gr_mat_randtest(A, state, ctx);
            status |= _gr_vec_randtest(cp, state, n + 1, ctx);
            status |= _gr_vec_randtest(cp2, state, n + 1, ctx);

            status |= _gr_mat_charpoly_berkowitz(cp, A, ctx);
            status |= charpoly_impl(cp2, A, ctx);

            /* Check that the output isn't just 0 */
            if (status == GR_SUCCESS && _gr_vec_equal(cp, cp2, n + 1, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                gr_ctx_println(ctx);
                flint_printf("A = "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("cp = "); _gr_vec_print(cp, n + 1, ctx); flint_printf("\n");
                flint_printf("cp2 = "); _gr_vec_print(cp2, n + 1, ctx); flint_printf("\n");
                flint_abort();
            }

            gr_mat_clear(A, ctx);
            gr_heap_clear_vec(cp, n + 1, ctx);
            gr_heap_clear_vec(cp2, n + 1, ctx);
        }

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}

POP_OPTIONS
