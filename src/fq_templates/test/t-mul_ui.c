/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Andres Goens
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
#include "long_extras.h"
#include "fmpz_poly.h"

TEST_TEMPLATE_FUNCTION_START(T, mul_ui, state)
{
    int i, result;

    /* Check aliasing of a, b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        ulong x;
        TEMPLATE(T, t) a, b;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        x = z_randtest(state);
        TEMPLATE(T, mul_ui)(b, a, x, ctx);
        TEMPLATE(T, mul_ui)(a, a, x, ctx);

        result = (TEMPLATE(T, equal)(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("x = %wu\n",x);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* compare with direct multiplication */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        ulong x;
        TEMPLATE(T, t) a, c;
        fmpz_poly_t b;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(c, ctx);
        fmpz_poly_init(b);

        TEMPLATE(T, randtest)(a, state, ctx);
        x = z_randtest(state);
        TEMPLATE(T, mul_ui)(c, a, x, ctx);
        fmpz_poly_scalar_mul_ui(b,a,x);
        TEMPLATE(T, reduce)(b,ctx);

        result = (TEMPLATE(T, equal)(c, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("x = %wu\n",x);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(c, ctx);
        fmpz_poly_clear(b);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
