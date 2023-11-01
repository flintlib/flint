/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
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

TEST_TEMPLATE_FUNCTION_START(T, sqr, state)
{
    int i, result;

    /* Check aliasing: a = a * a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, t) a, c;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);

        TEMPLATE(T, sqr)(c, a, ctx);
        TEMPLATE(T, sqr)(a, a, ctx);

        result = (TEMPLATE(T, equal)(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check a^2 + a^2 = a(a + a) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c, d;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);
        TEMPLATE(T, init)(d, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);

        TEMPLATE(T, sqr)(b, a, ctx);
        TEMPLATE(T, add)(c, b, b, ctx);

        TEMPLATE(T, add)(d, a, a, ctx);
        TEMPLATE(T, mul)(d, a, d, ctx);

        result = (TEMPLATE(T, equal)(c, d, ctx));
        if (!result)
        {
            flint_printf("FAIL (a^2 + a^2 == a(a + a)):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            flint_printf("d = "), TEMPLATE(T, print_pretty)(d, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);
        TEMPLATE(T, clear)(d, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
