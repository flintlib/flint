/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("mul_si....");
    fflush(stdout);

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
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_nmod_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
	    flint_printf("x = %wd\n",x);
            abort();
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
            nmod_poly_scalar_mul_nmod(b, a, -x);
            nmod_poly_neg(b, b);
        }
        else
        {
            nmod_poly_scalar_mul_nmod(b, a, x);
        }

        result = (fq_nmod_equal(c, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), fq_nmod_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
	    flint_printf("x = %wd\n",x);
            abort();
        }

        fq_nmod_clear(a, ctx);
        fq_nmod_clear(c, ctx);
        nmod_poly_clear(b);
        fq_nmod_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
