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

TEST_TEMPLATE_FUNCTION_START(T, is_invertible_f, state)
{
    slong ix;
    int result;

    /*
        Compare with the gcdinv function.

        N.B.  I checked by hand that this test shows both outcomes,
        i.e. trivial and non-trivial factors, sufficiently frequently.
     */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, ainv, f, g;

#if defined(FQ_ZECH_H)
        TEMPLATE(T, ctx_init_randtest_reducible)(ctx, state, 1);
#else
        TEMPLATE(T, ctx_init_randtest_reducible)(ctx, state, 0);
#endif

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(f, ctx);
        TEMPLATE(T, init)(g, ctx);
        TEMPLATE(T, init)(ainv, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, gcdinv)(f, ainv, a, ctx);

        result = (TEMPLATE(T, is_one)(f, ctx) ==
                  TEMPLATE(T, is_invertible_f)(g, a, ctx));
        result = result && TEMPLATE(T, equal)(f, g, ctx);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx);
            flint_printf("\n\n");
            flint_printf("f = "), TEMPLATE(T, print_pretty)(f, ctx);
            flint_printf("\n\n");
            flint_printf("g = "), TEMPLATE(T, print_pretty)(g, ctx);
            flint_printf("\n\n");
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(ainv, ctx);
        TEMPLATE(T, clear)(f, ctx);
        TEMPLATE(T, clear)(g, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
