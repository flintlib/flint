/*
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

TEST_TEMPLATE_FUNCTION_START(T, poly_get_str_pretty, state)
{
    int i, len;
    char *str;
    TEMPLATE(T, poly_t) a;
    TEMPLATE(T, ctx_t) ctx;

    TEMPLATE(T, ctx_randtest) (ctx, state);

    TEMPLATE(T, poly_init) (a, ctx);
    for (len = 0; len < 100; len++)
        for (i = 0; i < 1 * flint_test_multiplier(); i++)
        {
            TEMPLATE(T, poly_randtest) (a, state, len, ctx);
            str = TEMPLATE(T, poly_get_str_pretty) (a, "x", ctx);
            /* flint_printf("\n\n"); */
            /* TEMPLATE(T, poly_print_pretty)(a, "x", ctx); */
            /* flint_printf("\n%s\n", str); */
            flint_free(str);
        }

    TEMPLATE(T, poly_clear) (a, ctx);
    TEMPLATE(T, ctx_clear) (ctx);

    TEST_FUNCTION_END(state);
}
#endif
