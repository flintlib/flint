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

void _gr_poly_test_xgcd(gr_method_poly_xgcd_op xgcd_impl,
    flint_rand_t state, slong iters, slong maxn, gr_ctx_t ctx)
{
    slong iter;
    gr_ctx_ptr given_ctx = ctx;

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t my_ctx;
        gr_ctx_struct * ctx;

        if (given_ctx == NULL)
        {
            gr_ctx_init_random(my_ctx, state);
            ctx = my_ctx;
        }
        else
            ctx = given_ctx;

        {
            gr_poly_t a, b, d, g, s, t, v, w;
            int status = GR_SUCCESS;
            int aliasing;

            gr_poly_init(a, ctx);
            gr_poly_init(b, ctx);
            gr_poly_init(d, ctx);
            gr_poly_init(g, ctx);
            gr_poly_init(s, ctx);
            gr_poly_init(t, ctx);
            gr_poly_init(v, ctx);
            gr_poly_init(w, ctx);

            status |= gr_poly_randtest(a, state, 1 + n_randint(state, maxn), ctx);
            status |= gr_poly_randtest(b, state, 1 + n_randint(state, maxn), ctx);

            /* common factor */
            if (n_randint(state, 2))
            {
                status |= gr_poly_randtest(t, state, 1 + n_randint(state, maxn), ctx);
                status |= gr_poly_mul(a, a, t, ctx);
                status |= gr_poly_mul(b, b, t, ctx);
            }

            status |= gr_poly_randtest(g, state, 3, ctx);
            status |= gr_poly_randtest(s, state, 3, ctx);
            status |= gr_poly_randtest(t, state, 3, ctx);

            aliasing = n_randint(state, 8);

            switch (aliasing)
            {
                case 0:
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, b, ctx);
                    break;
                case 1:
                    status |= gr_poly_set(g, a, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, g, b, ctx);
                    break;
                case 2:
                    status |= gr_poly_set(s, a, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, s, b, ctx);
                    break;
                case 3:
                    status |= gr_poly_set(t, a, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, t, b, ctx);
                    break;
                case 4:
                    status |= gr_poly_set(g, b, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, g, ctx);
                    break;
                case 5:
                    status |= gr_poly_set(s, b, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, s, ctx);
                    break;
                case 6:
                    status |= gr_poly_set(t, b, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, t, ctx);
                    break;
                case 7:
                    status |= gr_poly_set(b, a, ctx);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, a, ctx);
                    break;
                default:
                    break;
            }

            status |= gr_poly_gcd_euclidean(d, a, b, ctx);

            status |= gr_poly_mul(v, s, a, ctx);
            status |= gr_poly_mul(w, t, b, ctx);
            status |= gr_poly_add(w, v, w, ctx);

            if (status == GR_SUCCESS && (gr_poly_equal(d, g, ctx) == T_FALSE || gr_poly_equal(g, w, ctx) == T_FALSE))
            {
                flint_printf("FAIL:\n");
                gr_poly_print(a, ctx), flint_printf("\n\n");
                gr_poly_print(b, ctx), flint_printf("\n\n");
                gr_poly_print(d, ctx), flint_printf("\n\n");
                gr_poly_print(g, ctx), flint_printf("\n\n");
                gr_poly_print(s, ctx), flint_printf("\n\n");
                gr_poly_print(t, ctx), flint_printf("\n\n");
                gr_poly_print(v, ctx), flint_printf("\n\n");
                gr_poly_print(w, ctx), flint_printf("\n\n");
                flint_abort();
            }

            if ((ctx->which_ring == GR_CTX_FMPQ || (ctx->which_ring == GR_CTX_NMOD8 && gr_ctx_is_field(ctx) == T_TRUE)) && status != GR_SUCCESS)
            {
                flint_printf("FAIL: did not succeed over Q or Z/pZ\n\n");
                gr_poly_print(a, ctx), flint_printf("\n\n");
                gr_poly_print(b, ctx), flint_printf("\n\n");
                flint_abort();
            }

            gr_poly_clear(a, ctx);
            gr_poly_clear(b, ctx);
            gr_poly_clear(d, ctx);
            gr_poly_clear(g, ctx);
            gr_poly_clear(s, ctx);
            gr_poly_clear(t, ctx);
            gr_poly_clear(v, ctx);
            gr_poly_clear(w, ctx);
        }

        if (given_ctx == NULL)
            gr_ctx_clear(ctx);
    }
}
