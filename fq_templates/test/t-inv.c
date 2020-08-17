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

    flint_printf("inv... ");
    fflush(stdout);

    /* Check aliasing: a = ~a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest_not_zero)(a, state, ctx);
        TEMPLATE(T, set)(b, a, ctx);

        TEMPLATE(T, inv)(c, b, ctx);
        TEMPLATE(T, inv)(b, b, ctx);

        result = (TEMPLATE(T, equal)(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            TEMPLATE(T, ctx_print)(ctx);
            abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check a * ~a == 1 for units*/
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest_not_zero)(a, state, ctx);

        TEMPLATE(T, inv)(b, a, ctx);
        TEMPLATE(T, mul)(c, a, b, ctx);

        result = (TEMPLATE(T, is_one)(c, ctx));
        if (!result)
        {
            flint_printf("FAIL (a * (~a) == 1):\n\n");
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

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}



#endif
