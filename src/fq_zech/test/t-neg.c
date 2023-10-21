/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_zech.h"

TEST_FUNCTION_START(fq_zech_neg, state)
{
    int i, j, result;
    fq_zech_ctx_t ctx;

    for (j = 0; j < 10; j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        /* Check aliasing: a = -a */
        for (i = 0; i < 200; i++)
        {
            fq_zech_t a, b, c;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fq_zech_init(c, ctx);

            fq_zech_randtest(a, state, ctx);
            fq_zech_set(b, a, ctx);

            fq_zech_neg(c, b, ctx);
            fq_zech_neg(b, b, ctx);

            result = (fq_zech_equal(b, c, ctx));
            if (!result)
            {
                flint_printf("FAIL a = -a:\n\n");
                fq_zech_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), fq_zech_print_pretty(c, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            fq_zech_clear(c, ctx);
        }

        /* Check a - b == a + (-b) */
        for (i = 0; i < 2000; i++)
        {
            fq_zech_t a, b, c1, c2;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fq_zech_init(c1, ctx);
            fq_zech_init(c2, ctx);

            fq_zech_randtest(a, state, ctx);
            fq_zech_randtest(b, state, ctx);

            fq_zech_sub(c1, a, b, ctx);
            fq_zech_neg(c2, b, ctx);
            fq_zech_add(c2, a, c2, ctx);

            result = (fq_zech_equal(c1, c2, ctx));
            if (!result)
            {
                flint_printf("FAIL  a - b == a + (-b):\n\n");
                fq_zech_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c1 = "), fq_zech_print_pretty(c1, ctx), flint_printf("\n");
                flint_printf("c2 = "), fq_zech_print_pretty(c2, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            fq_zech_clear(c1, ctx);
            fq_zech_clear(c2, ctx);
        }

        fq_zech_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
