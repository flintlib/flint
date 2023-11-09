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

TEST_FUNCTION_START(fq_zech_inv, state)
{
    int j, i, result;
    fq_zech_ctx_t ctx;

    for (j = 0; j < 10; j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        /* Check aliasing: a = ~a */
        for (i = 0; i < 200; i++)
        {
            fq_zech_t a, b, c;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fq_zech_init(c, ctx);

            fq_zech_randtest_not_zero(a, state, ctx);
            fq_zech_set(b, a, ctx);

            fq_zech_inv(c, b, ctx);
            fq_zech_inv(b, b, ctx);

            result = (fq_zech_equal(b, c, ctx));
            if (!result)
            {
                flint_printf("FAIL (aliasing):\n\n");
                fq_zech_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), fq_zech_print_pretty(c, ctx), flint_printf("\n");
                fq_zech_ctx_print(ctx);
                fflush(stdout);
                flint_abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            fq_zech_clear(c, ctx);
        }

        /* Check a * ~a == 1 for units */
        for (i = 0; i < 200; i++)
        {
            fq_zech_t a, b, c;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fq_zech_init(c, ctx);

            fq_zech_randtest_not_zero(a, state, ctx);

            fq_zech_inv(b, a, ctx);
            fq_zech_mul(c, a, b, ctx);

            result = (fq_zech_is_one(c, ctx));
            if (!result)
            {
                flint_printf("FAIL (a * (~a) == 1):\n\n");
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

        fq_zech_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
