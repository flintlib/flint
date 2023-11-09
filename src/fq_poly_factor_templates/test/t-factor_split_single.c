/*
    Copyright (C) 2019 Edouard Rousseau

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_factor_split_single, state)
{
    int iter;

    /* Compute a random splitting polynomial then check factorization */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        int len;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, q, r;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (q, ctx);
        TEMPLATE(T, poly_init) (r, ctx);

        /* random splitting monic polynomial of degree len */

        TEMPLATE(T, poly_randtest_monic) (a, state, 2, ctx);

        while (TEMPLATE(T, poly_degree) (a, ctx) < len)
        {
            TEMPLATE(T, poly_randtest_monic) (b, state, 2, ctx);
            TEMPLATE(T, poly_mul)(a, a, b, ctx);
        }

        /* b should be a factor of a */
        TEMPLATE(T, poly_factor_split_single) (b, a, ctx);

        /* check that b divides a */
        TEMPLATE(T, poly_divrem) (q, r, a, b, ctx);
        if (!TEMPLATE(T, poly_is_zero) (r, ctx))
        {
            flint_printf("FAIL:\n");
            flint_printf
                ("Error: factor does not divide original polynomial\n");
            flint_printf("factor:\n");
            TEMPLATE(T, poly_print) (b, ctx);
            flint_printf("\n\n");
            flint_printf("polynomial:\n");
            TEMPLATE(T, poly_print) (a, ctx);
            flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (q, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
