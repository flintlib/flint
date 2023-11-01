/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "nmod_poly.h"
#include "fq_nmod.h"

TEST_FUNCTION_START(fq_nmod_mul_ui, state)
{
    int i, result;

    /* Check aliasing of a, b */
    for (i = 0; i < 2000; i++)
    {
        fq_nmod_ctx_t ctx;
        ulong x;
        fq_nmod_t a, b;

        fq_nmod_ctx_randtest(ctx, state);

        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);

        fq_nmod_randtest(a, state, ctx);
        x = z_randtest(state);
        fq_nmod_mul_ui(b, a, x, ctx);
        fq_nmod_mul_ui(a, a, x, ctx);

        result = (fq_nmod_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL 1:\n\n");
            flint_printf("a = "), fq_nmod_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("x = %wu\n",x);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);

        fq_nmod_ctx_clear(ctx);
    }

    /* compare with direct multiplication */
    for (i = 0; i < 2000; i++)
    {
        fq_nmod_ctx_t ctx;
        ulong x;
        fq_nmod_t a, c;
        nmod_poly_t b;

        fq_nmod_ctx_randtest(ctx, state);

        fq_nmod_init(a, ctx);
        fq_nmod_init(c, ctx);
        nmod_poly_init(b, ctx->mod.n);

        fq_nmod_randtest(a, state, ctx);
        x = n_randint(state, ctx->mod.n);
        fq_nmod_mul_ui(c, a, x, ctx);
        nmod_poly_scalar_mul_nmod(b,a,x);

        result = (fq_nmod_equal(c, b, ctx));
        if (!result)
        {
            flint_printf("FAIL 2:\n\n");
            flint_printf("a = "), fq_nmod_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), fq_nmod_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("x = %wu\n",x);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(c, ctx);
        nmod_poly_clear(b);
        fq_nmod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
