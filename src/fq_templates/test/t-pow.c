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

TEST_TEMPLATE_FUNCTION_START(T, pow, state)
{
    int i, result;

    /* Check aliasing: a = a^e */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b;
        fmpz_t e;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        fmpz_init(e);

        TEMPLATE(T, randtest)(a, state, ctx);
        fmpz_randtest_unsigned(e, state, 6);

        TEMPLATE(T, pow)(b, a, e, ctx);
        TEMPLATE(T, pow)(a, a, e, ctx);

        result = (TEMPLATE(T, equal)(a, b, ctx));
        if (!result)
        {
            flint_printf("FAIL (alias):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        fmpz_clear(e);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Compare with multiplication, for integral values */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, t) a, b, c;
        fmpz_t e, f;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);
        fmpz_init(f);
        fmpz_init(e);

        TEMPLATE(T, randtest)(a, state, ctx);
        fmpz_randtest_unsigned(e, state, 6);

        TEMPLATE(T, pow)(b, a, e, ctx);
        TEMPLATE(T, one)(c, ctx);
        for (fmpz_one(f); fmpz_cmp(f, e) <= 0; fmpz_add_ui(f, f, 1))
        {
            TEMPLATE(T, mul)(c, c, a, ctx);
        }

        result = (TEMPLATE(T, equal)(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (cmp with mul):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("e = "), fmpz_print(e), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);
        fmpz_clear(e);
        fmpz_clear(f);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
