/*
    Copyright (C) 2021, 2022 William Hart
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_sqrt, state)
{
    int i;

    /* Test aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b;
        int square1, square2;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, poly_init)(a, ctx);
        TEMPLATE(T, poly_init)(b, ctx);

        TEMPLATE(T, poly_randtest)(a, state, 1 + n_randint(state, 50), ctx);

        if (n_randint(state, 2))
            TEMPLATE(T, poly_mul)(a, a, a, ctx);

        square1 = TEMPLATE(T, poly_sqrt)(b, a, ctx);
        square2 = TEMPLATE(T, poly_sqrt)(a, a, ctx);

        if (square1 != square2 || (square1 && !TEMPLATE(T, poly_equal)(a, b, ctx)))
        {
            flint_printf("FAIL: aliasing:\n");
            flint_printf("square1 = %d, square2 = %d\n\n", square1, square2);
            flint_printf("a: "); TEMPLATE(T, poly_print)(a, ctx); flint_printf("\n\n");
            flint_printf("b: "); TEMPLATE(T, poly_print)(b, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, poly_clear)(b, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Test random squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, c;
        int square;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, poly_init)(a, ctx);
        TEMPLATE(T, poly_init)(b, ctx);
        TEMPLATE(T, poly_init)(c, ctx);

        TEMPLATE(T, poly_randtest)(a, state, 1 + n_randint(state, 50), ctx);
        TEMPLATE(T, poly_mul)(b, a, a, ctx);
        square = TEMPLATE(T, poly_sqrt)(c, b, ctx);

        if (!square)
        {
            flint_printf("FAIL: square reported nonsquare:\n");
            flint_printf("a: "); TEMPLATE(T, poly_print)(a, ctx); flint_printf("\n\n");
            flint_printf("b: "); TEMPLATE(T, poly_print)(b, ctx); flint_printf("\n\n");
            flint_printf("c: "); TEMPLATE(T, poly_print)(c, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_mul)(c, c, c, ctx);
        if (!TEMPLATE(T, poly_equal)(c, b, ctx))
        {
            flint_printf("FAIL: sqrt(b)^2 != b:\n");
            flint_printf("a: "); TEMPLATE(T, poly_print)(a, ctx); flint_printf("\n\n");
            flint_printf("b: "); TEMPLATE(T, poly_print)(b, ctx); flint_printf("\n\n");
            flint_printf("c: "); TEMPLATE(T, poly_print)(c, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, poly_clear)(b, ctx);
        TEMPLATE(T, poly_clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Test "almost" squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, c;
        slong j;
        int square;
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) t;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(t, ctx);

        TEMPLATE(T, poly_init)(a, ctx);
        TEMPLATE(T, poly_init)(b, ctx);
        TEMPLATE(T, poly_init)(c, ctx);

        TEMPLATE(T, poly_randtest_not_zero)(a, state, 1 + n_randint(state, 50), ctx);
        TEMPLATE(T, poly_mul)(b, a, a, ctx);

        j = n_randint(state, TEMPLATE(T, poly_length)(b, ctx));
        TEMPLATE(T, randtest)(t, state, ctx);

        TEMPLATE(T, set)(b->coeffs + j, t, ctx);
        _TEMPLATE(T, poly_normalise)(b, ctx);

        square = TEMPLATE(T, poly_sqrt)(c, b, ctx);

        if (square)
        {
            TEMPLATE(T, poly_mul)(c, c, c, ctx);
            if (!TEMPLATE(T, poly_equal)(c, b, ctx))
            {
                flint_printf("FAIL: sqrt(b)^2 != b:\n");
                flint_printf("a: "); TEMPLATE(T, poly_print)(a, ctx); flint_printf("\n\n");
                flint_printf("b: "); TEMPLATE(T, poly_print)(b, ctx); flint_printf("\n\n");
                flint_printf("c: "); TEMPLATE(T, poly_print)(c, ctx); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        TEMPLATE(T, poly_clear)(a, ctx);
        TEMPLATE(T, poly_clear)(b, ctx);
        TEMPLATE(T, poly_clear)(c, ctx);

        TEMPLATE(T, clear)(t, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
