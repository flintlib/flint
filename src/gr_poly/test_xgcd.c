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

    for (iter = 0; iter < iters; iter++)
    {
        gr_ctx_t ctx2;
        gr_ctx_struct * ctxptr;

        if (ctx == NULL)
        {
            gr_ctx_init_random(ctx2, state);
            ctxptr = ctx2;
        }
        else
            ctxptr = ctx;

        {
            gr_poly_t a, b, d, g, s, t, v, w;
            int status = GR_SUCCESS;
            int aliasing;

            gr_poly_init(a, ctxptr);
            gr_poly_init(b, ctxptr);
            gr_poly_init(d, ctxptr);
            gr_poly_init(g, ctxptr);
            gr_poly_init(s, ctxptr);
            gr_poly_init(t, ctxptr);
            gr_poly_init(v, ctxptr);
            gr_poly_init(w, ctxptr);

            status |= gr_poly_randtest(a, state, 1 + n_randint(state, maxn), ctxptr);
            status |= gr_poly_randtest(b, state, 1 + n_randint(state, maxn), ctxptr);

            /* common factor */
            if (n_randint(state, 2))
            {
                status |= gr_poly_randtest(t, state, 1 + n_randint(state, maxn), ctxptr);
                status |= gr_poly_mul(a, a, t, ctxptr);
                status |= gr_poly_mul(b, b, t, ctxptr);
            }

            status |= gr_poly_randtest(g, state, 3, ctxptr);
            status |= gr_poly_randtest(s, state, 3, ctxptr);
            status |= gr_poly_randtest(t, state, 3, ctxptr);

            aliasing = n_randint(state, 8);

            switch (aliasing)
            {
                case 0:
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, b, ctxptr);
                    break;
                case 1:
                    status |= gr_poly_set(g, a, ctxptr);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, g, b, ctxptr);
                    break;
                case 2:
                    status |= gr_poly_set(s, a, ctxptr);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, s, b, ctxptr);
                    break;
                case 3:
                    status |= gr_poly_set(t, a, ctxptr);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, t, b, ctxptr);
                    break;
                case 4:
                    status |= gr_poly_set(g, b, ctxptr);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, g, ctxptr);
                    break;
                case 5:
                    status |= gr_poly_set(s, b, ctxptr);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, s, ctxptr);
                    break;
                case 6:
                    status |= gr_poly_set(t, b, ctxptr);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, t, ctxptr);
                    break;
                case 7:
                    status |= gr_poly_set(b, a, ctxptr);
                    status |= gr_poly_xgcd_wrapper(xgcd_impl, g, s, t, a, a, ctxptr);
                    break;
                default:
                    break;
            }

            status |= gr_poly_gcd_euclidean(d, a, b, ctxptr);

            status |= gr_poly_mul(v, s, a, ctxptr);
            status |= gr_poly_mul(w, t, b, ctxptr);
            status |= gr_poly_add(w, v, w, ctxptr);

            if (status == GR_SUCCESS && (gr_poly_equal(d, g, ctxptr) == T_FALSE || gr_poly_equal(g, w, ctxptr) == T_FALSE))
            {
                flint_printf("FAIL:\n");
                gr_poly_print(a, ctxptr), flint_printf("\n\n");
                gr_poly_print(b, ctxptr), flint_printf("\n\n");
                gr_poly_print(d, ctxptr), flint_printf("\n\n");
                gr_poly_print(g, ctxptr), flint_printf("\n\n");
                gr_poly_print(s, ctxptr), flint_printf("\n\n");
                gr_poly_print(t, ctxptr), flint_printf("\n\n");
                gr_poly_print(v, ctxptr), flint_printf("\n\n");
                gr_poly_print(w, ctxptr), flint_printf("\n\n");
                flint_abort();
            }

            if ((ctxptr->which_ring == GR_CTX_FMPQ || (ctxptr->which_ring == GR_CTX_NMOD8 && gr_ctx_is_field(ctxptr) == T_TRUE)) && status != GR_SUCCESS)
            {
                flint_printf("FAIL: did not succeed over Q or Z/pZ\n\n");
                gr_poly_print(a, ctxptr), flint_printf("\n\n");
                gr_poly_print(b, ctxptr), flint_printf("\n\n");
                flint_abort();
            }

            gr_poly_clear(a, ctxptr);
            gr_poly_clear(b, ctxptr);
            gr_poly_clear(d, ctxptr);
            gr_poly_clear(g, ctxptr);
            gr_poly_clear(s, ctxptr);
            gr_poly_clear(t, ctxptr);
            gr_poly_clear(v, ctxptr);
            gr_poly_clear(w, ctxptr);
        }

        if (ctx == NULL)
            gr_ctx_clear(ctxptr);
    }
}
