/*
    Copyright (C) 2018 Luca De Feo

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

TEST_TEMPLATE_FUNCTION_START(T, multiplicative_order, state)
{
    int i, result;

    /* Test that the computed multiplicative order is a multiple of the real one */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, tmp;
        fmpz_t ord;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(tmp, ctx);
        TEMPLATE(T, randtest)(a, state, ctx);

        fmpz_init(ord);
        result = TEMPLATE(T, multiplicative_order)(ord, a, ctx);
        TEMPLATE(T, pow)(tmp, a, ord, ctx);

        if (result && !TEMPLATE(T, is_one)(tmp, ctx))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("ord = "), fmpz_print(ord), flint_printf("\n");
            TEMPLATE(T, ctx_print)(ctx);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(ord);
        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(tmp, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Test that the computed multiplicative order is coherent with powering by p-1 */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) X, a;
        fmpz_t ord, size, pm1;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(X, ctx);
        TEMPLATE(T, init)(a, ctx);
        fmpz_init(ord);
        fmpz_init(size);
        fmpz_init(pm1);

        TEMPLATE(T, gen)(X, ctx);
        if (TEMPLATE(T, is_primitive)(X, ctx))
        {
            fmpz_sub_ui(pm1, TEMPLATE(T, ctx_prime)(ctx), 1);
            TEMPLATE(T, pow)(a, X, pm1, ctx);
            result = TEMPLATE(T, multiplicative_order)(ord, a, ctx);
            fmpz_mul(ord, ord, pm1);

            TEMPLATE(T, ctx_order)(size, ctx);
            fmpz_sub(size, size, ord);

            if (result && !fmpz_is_one(size))
            {
                flint_printf("FAIL:\n\n");
                flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                flint_printf("ord = "), fmpz_print(ord), flint_printf("\n");
                TEMPLATE(T, ctx_print)(ctx);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(pm1);
        fmpz_clear(ord);
        fmpz_clear(size);
        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(X, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
