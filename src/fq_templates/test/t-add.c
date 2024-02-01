/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, add, state)
{
    slong ix;
    int result;

#if defined(FQ_ZECH_H)
    for (ix = 0; ix < 10 * flint_test_multiplier(); ix++)
#else
    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
#endif
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c, x, y;
        int type;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 0);

#if defined(FQ_ZECH_H)
        for (slong jx = 0; jx < 100; jx++)
        {
#endif
        type = n_randint(state, 5);

        TEMPLATE(T, init)(a, ctx);
        if (type != 2)
            TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);
        if (type >= 3)
            TEMPLATE(T, init)(x, ctx);
        if (type == 4)
            TEMPLATE(T, init)(y, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        if (type != 2)
            TEMPLATE(T, randtest)(b, state, ctx);

        if (type == 0)
        {
            /* Check aliasing: a = a + b */
            TEMPLATE(T, add)(c, a, b, ctx);
            TEMPLATE(T, add)(a, a, b, ctx);
        }
        else if (type == 1)
        {
            /* Check aliasing: b = a + b */
            TEMPLATE(T, add)(c, b, a, ctx);
            TEMPLATE(T, add)(a, b, a, ctx);
        }
        else if (type == 2)
        {
            /* Check aliasing: a = a + a */
            TEMPLATE(T, add)(c, a, a, ctx);
            TEMPLATE(T, add)(a, a, a, ctx);
        }
        else if (type == 3)
        {
            /* Check that a + b == b + a */
            TEMPLATE(T, add)(c, a, b, ctx);
            TEMPLATE(T, add)(x, b, a, ctx);
            FLINT_SWAP(TEMPLATE(T, struct), *a, *x);
        }
        else
        {
            /* Check that (a + b) + c == a + (b + c) */
            TEMPLATE(T, randtest)(c, state, ctx);

            TEMPLATE(T, add)(x, a, b, ctx);
            TEMPLATE(T, add)(x, x, c, ctx);
            TEMPLATE(T, add)(y, b, c, ctx);
            TEMPLATE(T, add)(y, a, y, ctx);
            FLINT_SWAP(TEMPLATE(T, struct), *a, *x);
            FLINT_SWAP(TEMPLATE(T, struct), *c, *y);
        }

        result = (TEMPLATE(T, equal)(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("type = %d\n", type);
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            if (type != 2)
                flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            if (type >= 3)
                flint_printf("x = "), TEMPLATE(T, print_pretty)(x, ctx), flint_printf("\n");
            if (type == 4)
                flint_printf("y = "), TEMPLATE(T, print_pretty)(y, ctx), flint_printf("\n");
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        if (type != 2)
            TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);
        if (type >= 3)
            TEMPLATE(T, clear)(x, ctx);
        if (type == 4)
            TEMPLATE(T, clear)(y, ctx);

#if defined(FQ_ZECH_H)
        }
#endif
        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
