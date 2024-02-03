/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_zech.h"

TEST_FUNCTION_START(fq_zech_inv, state)
{
    slong ix, jx;
    int result;

    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fq_zech_ctx_t ctx;

        fq_zech_ctx_init_randtest(ctx, state, 1);

        for (jx = 0; jx < 10; jx++)
        {
            fq_zech_t a, b, c;
            int type;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fq_zech_init(c, ctx);

            fq_zech_randtest_not_zero(a, state, ctx);

            type = n_randint(state, 2);

            if (type == 0)
            {
                /* Check aliasing: a = ~a */
                fq_zech_set(b, a, ctx);
                fq_zech_inv(c, b, ctx);
                fq_zech_inv(b, b, ctx);

                result = (fq_zech_equal(b, c, ctx));
            }
            else
            {
                /* Check a * ~a == 1 for units */
                fq_zech_inv(b, a, ctx);
                fq_zech_mul(c, a, b, ctx);

                result = (fq_zech_is_one(c, ctx));
            }

            if (!result)
            {
                flint_printf("FAIL:\n\n");
                fq_zech_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), fq_zech_print_pretty(c, ctx), flint_printf("\n");
                fq_zech_ctx_print(ctx);
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
