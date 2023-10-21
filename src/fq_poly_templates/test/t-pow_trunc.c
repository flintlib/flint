/*
    Copyright (C) 2010, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_pow_trunc, state)
{
    int i, result;

    /* Check powering against naive method */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, c;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) d;
        slong e, trunc;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, init) (d, ctx);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 30), ctx);
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        TEMPLATE(T, poly_pow_trunc) (b, a, e, trunc, ctx);

        TEMPLATE(T, poly_pow) (c, a, e, ctx);
        TEMPLATE(T, poly_truncate) (c, trunc, ctx);

        TEMPLATE(T, poly_get_coeff) (d, c, 0, ctx);
        result = (TEMPLATE(T, poly_equal) (b, c, ctx)
                || (a->length == 0 && e == 0 && c->length == 1 && TEMPLATE(T, is_one) (d, ctx)));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, exp = %wd, trunc = %wd\n",
                    a->length, e, trunc);
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear) (d, ctx);
        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, c;
        TEMPLATE(T, ctx_t) ctx;
        slong e, trunc;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 30), ctx);
        e = n_randint(state, 20);
        trunc = n_randint(state, 30);

        TEMPLATE(T, poly_pow_trunc) (b, a, e, trunc, ctx);

        TEMPLATE(T, poly_set) (c, a, ctx);
        TEMPLATE(T, poly_pow_trunc) (c, c, e, trunc, ctx);

        result = (TEMPLATE(T, poly_equal) (b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, exp = %wd, trunc = %wd\n",
                    a->length, e, trunc);
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print) (c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
