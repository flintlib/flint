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

TEST_TEMPLATE_FUNCTION_START(T, is_square, state)
{
    int i, result;

    /* Check is_square(a^2) == 1 */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);

        TEMPLATE(T, sqr)(b, a, ctx);

        result = (TEMPLATE(T, is_square)(b, ctx));
        if (!result)
        {
            flint_printf("FAIL (is_square(a^2)):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check non-squares */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, z;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) != 0)
        {
            TEMPLATE(T, init)(a, ctx);
            TEMPLATE(T, init)(b, ctx);
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
                fflush(stdout);
                flint_abort();
            }

            TEMPLATE(T, clear)(a, ctx);
            TEMPLATE(T, clear)(b, ctx);
            TEMPLATE(T, clear)(z, ctx);

            TEMPLATE(T, ctx_clear)(ctx);
        }
    }

    TEST_FUNCTION_END(state);
}
#endif
