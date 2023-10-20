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

TEST_TEMPLATE_FUNCTION_START(T, is_primitive, state)
{
    int i;

    /* Test that is_primitive gives consistent answers with multiplicative_order */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a;
        fmpz_t ord, field_ord;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, randtest)(a, state, ctx);

        fmpz_init(ord);
        fmpz_init(field_ord);
        TEMPLATE(T, multiplicative_order)(ord, a, ctx);
        TEMPLATE(T, ctx_order)(field_ord, ctx);
        fmpz_sub(field_ord, field_ord, ord);

        if (TEMPLATE(T, is_primitive)(a, ctx) != fmpz_is_one(field_ord))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            TEMPLATE(T, ctx_print)(ctx);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(ord);
        fmpz_clear(field_ord);
        TEMPLATE(T, clear)(a, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}

#endif
