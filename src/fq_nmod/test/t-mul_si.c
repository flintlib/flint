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
#include "fmpz.h"
#include "fq_nmod.h"

TEST_FUNCTION_START(fq_nmod_mul_si, state)
{
    int i, result;

    {
        fq_nmod_t rop;
        fq_nmod_ctx_t ctx;
        fmpz_t p;
        fmpz_t f;

        fmpz_init(f);
        fmpz_init_set_si(p, 23);
        fq_nmod_ctx_init(ctx, p, 1, "x");

        fq_nmod_init(rop, ctx);

        for (i = 8; i <= 10; i++)
        {
            fq_nmod_set_si(rop, 1, ctx);
            fq_nmod_mul_si(rop, rop, 23 + i, ctx);
            fq_nmod_get_fmpz(f, rop, ctx);
            if (!fmpz_equal_si(f, i))
            {
                flint_printf("FAIL: f should have been %wd\n", i);
                flint_printf("f = "), fmpz_print(f), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_clear(rop, ctx);

        fq_nmod_ctx_clear(ctx);
        fmpz_clear(p);
        fmpz_clear(f);
    }

    /* Check aliasing of a, b */
    for (i = 0; i < 1000; i++)
    {
        fq_nmod_ctx_t ctx;
        slong x;
        fq_nmod_t a, b;

        fq_nmod_ctx_randtest(ctx, state);

        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);

        fq_nmod_randtest(a, state, ctx);
        x = z_randtest(state);
        fq_nmod_mul_si(b, a, x, ctx);
        fq_nmod_mul_si(a, a, x, ctx);

        result = (fq_nmod_equal(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL: check aliasing\n");
            flint_printf("a = "), fq_nmod_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
	        flint_printf("x = %wd\n",x);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);

        fq_nmod_ctx_clear(ctx);
    }

    /* compare with direct multiplication */
    for (i = 0; i < 1000; i++)
    {
        fq_nmod_ctx_t ctx;
        slong x;
        fq_nmod_t a, c;
	    nmod_poly_t b;

        fq_nmod_ctx_randtest(ctx, state);

        fq_nmod_init(a, ctx);
        fq_nmod_init(c, ctx);
	    nmod_poly_init(b, ctx->mod.n);

        fq_nmod_randtest(a, state, ctx);
        x = z_randtest(state);
        fq_nmod_mul_si(c, a, x, ctx);
        if (x < 0)
        {
            nmod_poly_scalar_mul_nmod(b, a, nmod_set_ui(-x, ctx->mod));
            nmod_poly_neg(b, b);
        }
        else
        {
            nmod_poly_scalar_mul_nmod(b, a, nmod_set_ui(x, ctx->mod));
        }

        result = (fq_nmod_equal(c, b, ctx));
        if (!result)
        {
            flint_printf("FAIL: check direct\n");
            flint_printf("a = "), fq_nmod_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
	        flint_printf("x = %wd\n",x);
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
