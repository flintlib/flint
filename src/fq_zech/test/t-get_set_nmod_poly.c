/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly.h"
#include "fq_zech.h"

TEST_FUNCTION_START(fq_zech_get_set_nmod_poly, state)
{
    slong i, j;
    fq_zech_ctx_t ctx;

    for (j = 0; j < 10*flint_test_multiplier(); j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        for (i = 0; i < 100; i++)
        {
            fq_zech_t a, b;
            nmod_poly_t c, t;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            nmod_poly_init(c, 1); /* modulus don't care */
            nmod_poly_init_mod(t, ctx->fq_nmod_ctx->modulus->mod);

            fq_zech_randtest(a, state, ctx);
            fq_zech_get_nmod_poly(c, a, ctx);
            nmod_poly_randtest(t, state, 20);
            nmod_poly_mul(t, t, ctx->fq_nmod_ctx->modulus);
            nmod_poly_add(c, c, t);
            fq_zech_set_nmod_poly(b, c, ctx);

            if (!fq_zech_equal(a, b, ctx))
            {
                flint_printf("FAIL:n\n");
                fq_zech_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), nmod_poly_print_pretty(c, "x"), flint_printf("\n");
                flint_printf("table = %wd\n", ctx->eval_table[a->value]);
                fflush(stdout);
                flint_abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            nmod_poly_clear(c);
            nmod_poly_clear(t);
        }

        fq_zech_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
