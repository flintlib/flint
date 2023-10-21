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

TEST_TEMPLATE_FUNCTION_START(T, poly_randtest_irreducible, state)
{
    int iter;

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) poly;
        slong length;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (poly, ctx);

        length = n_randint(state, 20) + 2;
        TEMPLATE(T, poly_randtest_irreducible) (poly, state, length, ctx);

        if (!TEMPLATE(T, poly_is_irreducible) (poly, ctx))
        {
            flint_printf("Error: reducible polynomial created!\n");
            flint_printf("poly:\n");
            TEMPLATE(T, poly_print_pretty) (poly, "x", ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (poly, ctx);

        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
