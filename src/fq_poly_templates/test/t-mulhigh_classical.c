/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_mulhigh_classical, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c;
        slong j, start;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest) (c, state, n_randint(state, 50), ctx);
        start = n_randint(state, 50);

        TEMPLATE(T, poly_mulhigh_classical) (a, b, c, start, ctx);
        TEMPLATE(T, poly_mulhigh_classical) (b, b, c, start, ctx);

        for (j = 0; j < start; j++)
        {
            if (j < a->length)
                TEMPLATE(T, zero) (a->coeffs + j, ctx);
            if (j < b->length)
                TEMPLATE(T, zero) (b->coeffs + j, ctx);
        }
        _TEMPLATE(T, poly_normalise) (a, ctx);
        _TEMPLATE(T, poly_normalise) (b, ctx);

        result = (TEMPLATE(T, poly_equal) (a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, poly_print_pretty) (a, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (b, "x", ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c;
        slong j, start;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest) (c, state, n_randint(state, 50), ctx);
        start = n_randint(state, 50);

        TEMPLATE(T, poly_mulhigh_classical) (a, b, c, start, ctx);
        TEMPLATE(T, poly_mulhigh_classical) (c, b, c, start, ctx);

        for (j = 0; j < start; j++)
        {
            if (j < a->length)
                TEMPLATE(T, zero) (a->coeffs + j, ctx);
            if (j < c->length)
                TEMPLATE(T, zero) (c->coeffs + j, ctx);
        }

        _TEMPLATE(T, poly_normalise) (a, ctx);
        _TEMPLATE(T, poly_normalise) (c, ctx);

        result = (TEMPLATE(T, poly_equal) (a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, poly_print_pretty) (a, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (c, "x", ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Compare with mul_basecase */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, d;
        slong j, start;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest) (c, state, n_randint(state, 50), ctx);
        start = n_randint(state, 50);

        TEMPLATE(T, poly_mul_classical) (a, b, c, ctx);
        TEMPLATE(T, poly_mulhigh_classical) (d, b, c, start, ctx);

        for (j = 0; j < start; j++)
        {
            if (j < a->length)
                TEMPLATE(T, zero) (a->coeffs + j, ctx);
            if (j < d->length)
                TEMPLATE(T, zero) (d->coeffs + j, ctx);
        }

        _TEMPLATE(T, poly_normalise) (a, ctx);
        _TEMPLATE(T, poly_normalise) (d, ctx);

        result = (TEMPLATE(T, poly_equal) (a, d, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            TEMPLATE(T, poly_print_pretty) (b, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (c, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (a, "x", ctx), flint_printf("\n\n");
            TEMPLATE(T, poly_print_pretty) (d, "x", ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
