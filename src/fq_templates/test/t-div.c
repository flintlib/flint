/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, div, state)
{
    int i, result;

    /* Check aliasing: a = a * b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c, d;

#if defined(FQ_ZECH_H)
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 1);
#else
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 0);
#endif

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);
        TEMPLATE(T, init)(d, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, randtest_not_zero)(b, state, ctx);

        TEMPLATE(T, div)(c, a, b, ctx);
        TEMPLATE(T, mul)(d, c, b, ctx);

        result = (TEMPLATE(T, equal)(a, d, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
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
