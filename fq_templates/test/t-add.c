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

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("add... ");
    fflush(stdout);

    /* Check aliasing: a = a + b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;

        TEMPLATE(T, ctx_randtest)(ctx, state);
        
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, randtest)(b, state, ctx);

        TEMPLATE(T, add)(c, a, b, ctx);
        TEMPLATE(T, add)(a, a, b, ctx);

        result = (TEMPLATE(T, equal)(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check aliasing: b = a + b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, randtest)(b, state, ctx);

        TEMPLATE(T, add)(c, a, b, ctx);
        TEMPLATE(T, add)(b, a, b, ctx);

        result = (TEMPLATE(T, equal)(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check aliasing: a = a + a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, c;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);

        TEMPLATE(T, add)(c, a, a, ctx);
        TEMPLATE(T, add)(a, a, a, ctx);

        result = (TEMPLATE(T, equal)(a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check that a + b == b + a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c1, c2;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c1, ctx);
        TEMPLATE(T, init)(c2, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, randtest)(b, state, ctx);

        TEMPLATE(T, add)(c1, a, b, ctx);
        TEMPLATE(T, add)(c2, b, a, ctx);

        result = (TEMPLATE(T, equal)(c1, c2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a  = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b  = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c1 = "), TEMPLATE(T, print_pretty)(c1, ctx), flint_printf("\n");
            flint_printf("c2 = "), TEMPLATE(T, print_pretty)(c2, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c1, ctx);
        TEMPLATE(T, clear)(c2, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check that (a + b) + c == a + (b + c) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c, lhs, rhs;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);
        TEMPLATE(T, init)(lhs, ctx);
        TEMPLATE(T, init)(rhs, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, randtest)(b, state, ctx);
        TEMPLATE(T, randtest)(c, state, ctx);

        TEMPLATE(T, add)(lhs, a, b, ctx);
        TEMPLATE(T, add)(lhs, lhs, c, ctx);
        TEMPLATE(T, add)(rhs, b, c, ctx);
        TEMPLATE(T, add)(rhs, a, rhs, ctx);

        result = (TEMPLATE(T, equal)(lhs, rhs, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a   = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b   = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c   = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            flint_printf("lhs = "), TEMPLATE(T, print_pretty)(lhs, ctx), flint_printf("\n");
            flint_printf("rhs = "), TEMPLATE(T, print_pretty)(rhs, ctx), flint_printf("\n");
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);
        TEMPLATE(T, clear)(lhs, ctx);
        TEMPLATE(T, clear)(rhs, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}



#endif
