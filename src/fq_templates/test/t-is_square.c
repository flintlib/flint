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

TEST_TEMPLATE_FUNCTION_START(T, is_square, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 300 * flint_test_multiplier(); ix++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b;
        int type;

        type = n_randint(state, 2);

        if (type == 0)
#if defined(FQ_ZECH_H)
            TEMPLATE(T, ctx_init_randtest)(ctx, state, 1);
#else
            TEMPLATE(T, ctx_init_randtest)(ctx, state, 0);
#endif
        else
            do
#if defined(FQ_ZECH_H)
                TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
#else
                TEMPLATE(T, ctx_init_randtest)(ctx, state, 1);
#endif
#if defined(FQ_NMOD_H) || defined(FQ_ZECH_H)
            while (TEMPLATE(T, ctx_prime)(ctx) == 2);
#else
            while (fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) == 0);
#endif

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);

        if (type == 0)
        {
            /* Check is_square(a^2) == 1 */
            TEMPLATE(T, randtest)(a, state, ctx);
            TEMPLATE(T, sqr)(b, a, ctx);

            result = (TEMPLATE(T, is_square)(b, ctx));
            if (!result)
            {
                flint_printf("FAIL (is_square(a^2)):\n\n");
                flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
                flint_abort();
            }
        }
        else
        {
            /* Check non-squares */
            TEMPLATE(T, t) z;

            TEMPLATE(T, init)(z, ctx);

            while (TEMPLATE(T, is_square)(z, ctx))
                TEMPLATE(T, randtest)(z, state, ctx);

            while (TEMPLATE(T, is_zero)(a, ctx))
                TEMPLATE(T, randtest)(a, state, ctx);

            TEMPLATE(T, sqr)(b, a, ctx);
            TEMPLATE(T, mul)(b, b, z, ctx);

            result = (!TEMPLATE(T, is_square)(b, ctx));
            if (!result)
            {
                flint_printf("FAIL (is_square(a^2)):\n\n");
                flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                flint_printf("z = "), TEMPLATE(T, print_pretty)(z, ctx), flint_printf("\n");
                flint_abort();
            }

            TEMPLATE(T, clear)(z, ctx);
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
