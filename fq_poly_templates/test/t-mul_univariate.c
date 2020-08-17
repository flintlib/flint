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

    flint_printf("mul_univariate... ");
    fflush(stdout);

    /* Check aliasing: a = a * b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest) (b, state, len, ctx);

        TEMPLATE(T, poly_mul_univariate) (c, a, b, ctx);
        TEMPLATE(T, poly_mul_univariate) (a, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing: b = a * b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest) (b, state, len, ctx);

        TEMPLATE(T, poly_mul_univariate) (c, a, b, ctx);
        TEMPLATE(T, poly_mul_univariate) (b, a, b, ctx);

        result = (TEMPLATE(T, poly_equal) (b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check aliasing: a = a * a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, c;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (c, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);

        TEMPLATE(T, poly_mul_univariate) (c, a, a, ctx);
        TEMPLATE(T, poly_mul_univariate) (a, a, a, ctx);

        result = (TEMPLATE(T, poly_equal) (a, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check that a * b == b * a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a, b, c, e;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (e, ctx);

        TEMPLATE(T, poly_randtest) (a, state, len, ctx);
        TEMPLATE(T, poly_randtest) (b, state, len, ctx);

        TEMPLATE(T, poly_mul_univariate) (c, a, b, ctx);
        TEMPLATE(T, poly_mul_univariate) (e, b, a, ctx);

        result = (TEMPLATE(T, poly_equal) (e, c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, poly_print_pretty) (a, "X", ctx),
                flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, poly_print_pretty) (b, "X", ctx),
                flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, poly_print_pretty) (c, "X", ctx),
                flint_printf("\n");
            flint_printf("e = "), TEMPLATE(T, poly_print_pretty) (e, "X", ctx),
                flint_printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (e, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Check that (b*c)+(b*d) = b*(c+d) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        slong len;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, poly_t) a1, a2, b, c, d;

        len = n_randint(state, 15) + 1;
        TEMPLATE(T, ctx_randtest) (ctx, state);
        TEMPLATE(T, poly_init) (a1, ctx);
        TEMPLATE(T, poly_init) (a2, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);

        TEMPLATE(T, poly_randtest) (b, state, len, ctx);
        TEMPLATE(T, poly_randtest) (c, state, len, ctx);
        TEMPLATE(T, poly_randtest) (d, state, len, ctx);

        TEMPLATE(T, poly_mul_univariate) (a1, b, c, ctx);
        TEMPLATE(T, poly_mul_univariate) (a2, b, d, ctx);
        TEMPLATE(T, poly_add) (a1, a1, a2, ctx);

        TEMPLATE(T, poly_add) (c, c, d, ctx);
        TEMPLATE(T, poly_mul_univariate) (a2, b, c, ctx);

        result = (TEMPLATE(T, poly_equal) (a1, a2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a1 = "), TEMPLATE(T, poly_print_pretty) (a1, "X",
                                                                   ctx),
                flint_printf("\n");
            flint_printf("a2 = "), TEMPLATE(T, poly_print_pretty) (a2, "X",
                                                                   ctx),
                flint_printf("\n");
            flint_printf("b  = "), TEMPLATE(T, poly_print_pretty) (b, "X",
                                                                   ctx),
                flint_printf("\n");
            flint_printf("c  = "), TEMPLATE(T, poly_print_pretty) (c, "X",
                                                                   ctx),
                flint_printf("\n");
            flint_printf("d  = "), TEMPLATE(T, poly_print_pretty) (d, "X",
                                                                   ctx),
                flint_printf("\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a1, ctx);
        TEMPLATE(T, poly_clear) (a2, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}



#endif
