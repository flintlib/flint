/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"

TEST_TEMPLATE_FUNCTION_START(T, is_invertible, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 300 * flint_test_multiplier(); ix++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a;

#if defined(FQ_ZECH_H)
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);
#else
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 0);
#endif

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, randtest)(a, state, ctx);

        result = (TEMPLATE(T, is_invertible)(a, ctx) != TEMPLATE(T, is_zero)(a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            TEMPLATE(T, ctx_print)(ctx);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
