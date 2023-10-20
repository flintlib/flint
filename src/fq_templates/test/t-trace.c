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
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, trace, state)
{
    int i, result;

    /* Compare with sum of Galois conjugates */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;
        fmpz_t x, y;
        slong j;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);
        fmpz_init(x);
        fmpz_init(y);

        TEMPLATE(T, randtest)(a, state, ctx);

        TEMPLATE(T, trace)(x, a, ctx);

        TEMPLATE(T, zero)(b, ctx);
        for (j = 0; j < TEMPLATE(T, ctx_degree)(ctx); j++)
        {
            TEMPLATE(T, frobenius)(c, a, j, ctx);
            TEMPLATE(T, add)(b, b, c, ctx);
        }

        TEMPLATE(T, zero)(c, ctx);
        TEMPLATE(T, set_fmpz)(c, x, ctx);

        result = TEMPLATE(T, equal)(b, c, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            flint_printf("x = "), fmpz_print(x), flint_printf("\n");
            for (j = 0; j < TEMPLATE(T, ctx_degree)(ctx); j++)
            {
                TEMPLATE(T, frobenius)(c, a, j, ctx);
                flint_printf("sigma^%wd = ", j), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            }
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);
        fmpz_clear(x);
        fmpz_clear(y);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
