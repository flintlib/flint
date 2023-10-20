/*
    Copyright (C) 2020 William Hart

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

TEST_TEMPLATE_FUNCTION_START(T, sqrt, state)
{
    int j, i, result;
    TEMPLATE(T, ctx_t) ctx;

    for (j = 0; j < 10; j++)
    {
        TEMPLATE(T, ctx_randtest)(ctx, state);

        /* Check aliasing: a = a * a */
        for (i = 0; i < 200; i++)
        {
            TEMPLATE(T, t) a, c;

            TEMPLATE(T, init)(a, ctx);
            TEMPLATE(T, init)(c, ctx);

            TEMPLATE(T, randtest)(a, state, ctx);

            TEMPLATE(T, sqr)(c, a, ctx);

            TEMPLATE(T, sqrt)(a, c, ctx);
            TEMPLATE(T, sqrt)(c, c, ctx);

            result = (TEMPLATE(T, equal)(a, c, ctx));
            if (!result)
            {
                flint_printf("FAIL (aliasing):\n\n");
                flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            TEMPLATE(T, clear)(a, ctx);
            TEMPLATE(T, clear)(c, ctx);
        }

        /* Check sqrt(a^2) = a and that x*a^2 is not a square */
        for (i = 0; i < 200; i++)
        {
            int r;
            TEMPLATE(T, t) a, b, c, d, x;

            TEMPLATE(T, init)(a, ctx);
            TEMPLATE(T, init)(b, ctx);
            TEMPLATE(T, init)(c, ctx);
            TEMPLATE(T, init)(d, ctx);
            TEMPLATE(T, init)(x, ctx);

            TEMPLATE(T, randtest)(a, state, ctx);

            TEMPLATE(T, sqr)(b, a, ctx);

            r = TEMPLATE(T, sqrt)(c, b, ctx);
            TEMPLATE(T, sqr)(d, c, ctx);

            result = (r && TEMPLATE(T, equal)(d, b, ctx));
            if (!result)
            {
                flint_printf("FAIL (sqrt(a^2) == a):\n\n");
                flint_printf("r = %d\n", r);
                flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
                flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
                flint_printf("d = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            if (ctx->is_conway && fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) != 0 &&
                    !TEMPLATE(T, is_zero)(b, ctx))
            {
                TEMPLATE(T, gen)(x, ctx);
                TEMPLATE(T, mul)(b, b, x, ctx);

                r = TEMPLATE(T, sqrt)(c, b, ctx);

                result = !r; /* check b is not a square */
                if (!result)
                {
                    flint_printf("FAIL (a^2*x is nonsquare):\n\n");
                    flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                    flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            TEMPLATE(T, clear)(a, ctx);
            TEMPLATE(T, clear)(b, ctx);
            TEMPLATE(T, clear)(c, ctx);
            TEMPLATE(T, clear)(d, ctx);
            TEMPLATE(T, clear)(x, ctx);
        }

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
