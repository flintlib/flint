/*
    Copyright (C) 2009 William Hart
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

TEST_TEMPLATE_FUNCTION_START(T, poly_mulhigh, state)
{
    int i, result;

    /* Compare with left truncated product of a and b */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c;
        slong j, n;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        n = n_randint(state, 50);
        TEMPLATE(T, poly_randtest) (b, state, n, ctx);
        TEMPLATE(T, poly_randtest) (c, state, n, ctx);

        TEMPLATE(T, poly_mulhigh) (a, b, c, n, ctx);
        TEMPLATE(T, poly_mul) (b, b, c, ctx);
        for (j = 0; j < n; j++)
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

    TEST_FUNCTION_END(state);
}
#endif
