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

#include "fq_nmod.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    slong i, j;
    fq_nmod_ctx_t ctx;
    FLINT_TEST_INIT(state);
    
    flint_printf("get/set_nmod_poly... ");
    fflush(stdout);

    for (j = 0; j < 10*flint_test_multiplier(); j++)
    {
        fq_nmod_ctx_randtest(ctx, state);

        for (i = 0; i < 100; i++)
        {
            fq_nmod_t a, b;
            nmod_poly_t c, t;

            fq_nmod_init(a, ctx);
            fq_nmod_init(b, ctx);
            nmod_poly_init(c, 1); /* modulus don't care */
            nmod_poly_init_mod(t, ctx->modulus->mod);

            fq_nmod_randtest(a, state, ctx);
            fq_nmod_get_nmod_poly(c, a, ctx);
            nmod_poly_randtest(t, state, 20);
            nmod_poly_mul(t, t, ctx->modulus);
            nmod_poly_add(c, c, t);
            fq_nmod_set_nmod_poly(b, c, ctx);

            if (!fq_nmod_equal(a, b, ctx))
            {
                flint_printf("FAIL:n\n");
                fq_nmod_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("a = "), fq_nmod_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_nmod_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), nmod_poly_print_pretty(c, "x"), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_clear(a, ctx);
            fq_nmod_clear(b, ctx);
            nmod_poly_clear(c);
            nmod_poly_clear(t);
        }

        fq_nmod_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
