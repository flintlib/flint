/*
    Copyright (C) 2020 William Hart
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
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, sqrt, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 300 * flint_test_multiplier(); ix++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, c;
        int type;

        type = n_randint(state, 2);

#if defined(FQ_ZECH_H)
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
#else
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 1);
#endif

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(c, ctx);

        if (type == 0)
        {
            /* Check aliasing: a = a * a */
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
                flint_abort();
            }
        }
        else
        {
            /* Check sqrt(a^2) = a and that x*a^2 is not a square */
            TEMPLATE(T, t) b, d, x;

            TEMPLATE(T, init)(b, ctx);
            TEMPLATE(T, init)(d, ctx);
            TEMPLATE(T, init)(x, ctx);

            TEMPLATE(T, randtest)(a, state, ctx);

            TEMPLATE(T, sqr)(b, a, ctx);

            result = TEMPLATE(T, sqrt)(c, b, ctx);
            TEMPLATE(T, sqr)(d, c, ctx);

            result = (result && TEMPLATE(T, equal)(d, b, ctx));
            if (!result)
            {
                flint_printf("FAIL (sqrt(a^2) == a):\n\n");
                flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
                flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
                flint_printf("d = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
                flint_abort();
            }

#if defined(FQ_NMOD_H) || defined(FQ_ZECH_H)
            if (ctx->is_conway && TEMPLATE(T, ctx_prime)(ctx) != 2 &&
                    !TEMPLATE(T, is_zero)(b, ctx))
#else
            if (ctx->is_conway && fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) != 0 &&
                    !TEMPLATE(T, is_zero)(b, ctx))
#endif
            {
                TEMPLATE(T, gen)(x, ctx);
                TEMPLATE(T, mul)(b, b, x, ctx);

                result = !TEMPLATE(T, sqrt)(c, b, ctx); /* check b is not a square */
                if (!result)
                {
                    flint_printf("FAIL (a^2*x is nonsquare):\n\n");
                    flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                    flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
                    flint_abort();
                }
            }

            TEMPLATE(T, clear)(b, ctx);
            TEMPLATE(T, clear)(d, ctx);
            TEMPLATE(T, clear)(x, ctx);
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
