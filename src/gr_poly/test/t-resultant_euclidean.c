/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 William Hart
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

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_poly_resultant_euclidean, state)
{
    slong iter;

    /* Check res(f, g) == (-1)^(deg f deg g) res(g, f) */
    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t f, g;
        gr_ptr x, y;
        int status = GR_SUCCESS;
        slong n;

        gr_ctx_init_random(ctx, state);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        if (ctx->methods == _ca_methods)
            n = n_randint(state, 3);
        else
            n = n_randint(state, 6);

        status |= gr_poly_randtest(f, state, n, ctx);
        status |= gr_poly_randtest(g, state, n, ctx);
        status |= gr_poly_resultant_euclidean(x, f, g, ctx);
        status |= gr_poly_resultant_euclidean(y, g, f, ctx);

        if (((f->length - 1) * (g->length - 1)) % 2)
           status |= gr_neg(y, y, ctx);

        if (status == GR_SUCCESS && gr_equal(x, y, ctx) == T_FALSE)
        {
            flint_printf("FAIL (res(f, g) == (-1)^(deg f deg g) res(g, f)):\n");
            gr_ctx_println(ctx);
            gr_poly_print(f, ctx), flint_printf("\n\n");
            gr_poly_print(g, ctx), flint_printf("\n\n");
            printf("x = "); gr_println(x, ctx);
            printf("y = "); gr_println(y, ctx);
            fflush(stdout);
            flint_abort();
        }

        if ((ctx->which_ring == GR_CTX_FMPQ || (ctx->which_ring == GR_CTX_NMOD8 && gr_ctx_is_field(ctx) == T_TRUE)) && status != GR_SUCCESS)
        {
            flint_printf("FAIL: did not succeed over Q or Z/pZ\n\n");
            gr_ctx_println(ctx);
            gr_poly_print(f, ctx), flint_printf("\n\n");
            gr_poly_print(g, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);

        gr_ctx_clear(ctx);
    }

    /* Check res(f h, g) == res(f, g) res(h, g) */
    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_poly_t f, fh, h, g;
        gr_ptr x, y, z, yz;
        int status = GR_SUCCESS;
        slong n;

        gr_ctx_init_random(ctx, state);

        gr_poly_init(f, ctx);
        gr_poly_init(fh, ctx);
        gr_poly_init(g, ctx);
        gr_poly_init(h, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        z = gr_heap_init(ctx);
        yz = gr_heap_init(ctx);

        if (ctx->methods == _ca_methods)
            n = n_randint(state, 3);
        else
            n = n_randint(state, 6);

        status |= gr_poly_randtest(f, state, n, ctx);
        status |= gr_poly_randtest(g, state, n, ctx);
        status |= gr_poly_randtest(h, state, n, ctx);
        status |= gr_poly_resultant_euclidean(y, f, g, ctx);
        status |= gr_poly_resultant_euclidean(z, h, g, ctx);
        status |= gr_mul(yz, y, z, ctx);

        status |= gr_poly_mul(fh, f, h, ctx);
        status |= gr_poly_resultant_euclidean(x, fh, g, ctx);

        if (status == GR_SUCCESS && gr_ctx_is_field(ctx) == T_TRUE && gr_equal(x, yz, ctx) == T_FALSE)
        {
            flint_printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            gr_ctx_println(ctx);
            printf("f = "); gr_poly_print(f, ctx), flint_printf("\n\n");
            printf("g = "); gr_poly_print(g, ctx), flint_printf("\n\n");
            printf("h = "); gr_poly_print(h, ctx), flint_printf("\n\n");
            printf("fh = "); gr_poly_print(fh, ctx), flint_printf("\n\n");
            printf("x = "); gr_println(x, ctx);
            printf("y = "); gr_println(y, ctx);
            printf("z = "); gr_println(z, ctx);
            printf("yz = "); gr_println(yz, ctx);
            fflush(stdout);
            flint_abort();
        }

        if ((ctx->which_ring == GR_CTX_FMPQ || (ctx->which_ring == GR_CTX_NMOD8 && gr_ctx_is_field(ctx) == T_TRUE)) && status != GR_SUCCESS)
        {
            flint_printf("FAIL: did not succeed over Q or Z/pZ\n\n");
            gr_ctx_println(ctx);
            gr_poly_print(f, ctx), flint_printf("\n\n");
            gr_poly_print(g, ctx), flint_printf("\n\n");
            gr_poly_print(h, ctx), flint_printf("\n\n");
            gr_poly_print(fh, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(fh, ctx);
        gr_poly_clear(g, ctx);
        gr_poly_clear(h, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_heap_clear(z, ctx);
        gr_heap_clear(yz, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
