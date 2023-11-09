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

TEST_FUNCTION_START(fq_zech_sqr, state)
{
    int j, i, result;
    fq_zech_ctx_t ctx;

    for (j = 0; j < 10; j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        /* Check aliasing: a = a * a */
        for (i = 0; i < 200; i++)
        {
            fq_zech_t a, c;

            fq_zech_init(a, ctx);
            fq_zech_init(c, ctx);

            fq_zech_randtest(a, state, ctx);

            fq_zech_sqr(c, a, ctx);
            fq_zech_sqr(a, a, ctx);

            result = (fq_zech_equal(a, c, ctx));
            if (!result)
            {
                flint_printf("FAIL (aliasing):\n\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("c = "), fq_zech_print_pretty(c, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(c, ctx);
        }

        /* Check a^2 + a^2 = a(a + a) */
        for (i = 0; i < 200; i++)
        {
            fq_zech_t a, b, c, d;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fq_zech_init(c, ctx);
            fq_zech_init(d, ctx);

            fq_zech_randtest(a, state, ctx);

            fq_zech_sqr(b, a, ctx);
            fq_zech_add(c, b, b, ctx);

            fq_zech_add(d, a, a, ctx);
            fq_zech_mul(d, a, d, ctx);

            result = (fq_zech_equal(c, d, ctx));
            if (!result)
            {
                flint_printf("FAIL (a^2 + a^2 == a(a + a)):\n\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), fq_zech_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), fq_zech_print_pretty(d, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            fq_zech_clear(c, ctx);
            fq_zech_clear(d, ctx);

        }

        fq_zech_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
