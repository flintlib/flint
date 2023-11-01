/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_xgcd_hgcd, state)
{
    int i;

    /* Compare with result from GCD and check correctness */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        gr_ctx_t ctx;
        gr_poly_t a, b, d, g, s, t, v, w;
        int status = GR_SUCCESS;
        int aliasing;
        slong n, cutoff1, cutoff2;

        gr_ctx_init_random(ctx, state);

        cutoff1 = n_randint(state, 10);
        cutoff2 = n_randint(state, 10);

        if (gr_ctx_is_finite(ctx) == T_TRUE && n_randint(state, 2) == 0)
        {
            n = 100;
            cutoff1 = n_randint(state, 100);
            cutoff2 = n_randint(state, 100);
        }
        else if (ctx->which_ring == GR_CTX_CC_CA || ctx->which_ring == GR_CTX_RR_CA)
            n = 4;
        else
            n = 6;

        gr_poly_init(a, ctx);
        gr_poly_init(b, ctx);
        gr_poly_init(d, ctx);
        gr_poly_init(g, ctx);
        gr_poly_init(s, ctx);
        gr_poly_init(t, ctx);
        gr_poly_init(v, ctx);
        gr_poly_init(w, ctx);

        status |= gr_poly_randtest(a, state, n_randint(state, n), ctx);
        status |= gr_poly_randtest(b, state, n_randint(state, n), ctx);

        /* common factor */
        if (n_randint(state, 2))
        {
            status |= gr_poly_randtest(t, state, n_randint(state, n), ctx);
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
                status |= gr_poly_xgcd_hgcd(g, s, t, a, b, cutoff1, cutoff2, ctx);
                break;
            case 1:
                status |= gr_poly_set(g, a, ctx);
                status |= gr_poly_xgcd_hgcd(g, s, t, g, b, cutoff1, cutoff2, ctx);
                break;
            case 2:
                status |= gr_poly_set(s, a, ctx);
                status |= gr_poly_xgcd_hgcd(g, s, t, s, b, cutoff1, cutoff2, ctx);
                break;
            case 3:
                status |= gr_poly_set(t, a, ctx);
                status |= gr_poly_xgcd_hgcd(g, s, t, t, b, cutoff1, cutoff2, ctx);
                break;
            case 4:
                status |= gr_poly_set(g, b, ctx);
                status |= gr_poly_xgcd_hgcd(g, s, t, a, g, cutoff1, cutoff2, ctx);
                break;
            case 5:
                status |= gr_poly_set(s, b, ctx);
                status |= gr_poly_xgcd_hgcd(g, s, t, a, s, cutoff1, cutoff2, ctx);
                break;
            case 6:
                status |= gr_poly_set(t, b, ctx);
                status |= gr_poly_xgcd_hgcd(g, s, t, a, t, cutoff1, cutoff2, ctx);
                break;
            case 7:
                status |= gr_poly_set(b, a, ctx);
                status |= gr_poly_xgcd_hgcd(g, s, t, a, a, cutoff1, cutoff2, ctx);
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
            fflush(stdout);
            flint_abort();
        }

        if ((ctx->which_ring == GR_CTX_FMPQ || (ctx->which_ring == GR_CTX_NMOD8 && gr_ctx_is_field(ctx) == T_TRUE)) && status != GR_SUCCESS)
        {
            flint_printf("FAIL: did not succeed over Q or Z/pZ\n\n");
            gr_poly_print(a, ctx), flint_printf("\n\n");
            gr_poly_print(b, ctx), flint_printf("\n\n");
            fflush(stdout);
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
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
