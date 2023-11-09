/*
    Copyright (C) 2009, 2014 William Hart
    Copyright (C) 2011 Sebastian Pancratz
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

TEST_TEMPLATE_FUNCTION_START(T, poly_div_series, state)
{
    int i, result;

    /* Check A*B^{-1} * B is congruent A mod t^n */
    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, d;
        slong n = n_randint(state, 80) + 1;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (d, ctx);

        TEMPLATE(T, poly_randtest) (c, state, n_randint(state, 80) + 2, ctx);
        TEMPLATE(T, poly_randtest) (d, state, n_randint(state, 80) + 2, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 80) + 1, ctx);
        TEMPLATE(T, poly_randtest_not_zero) (b, state,
                                             n_randint(state, 80) + 1, ctx);
        TEMPLATE(T, randtest_not_zero) (b->coeffs + 0, state, ctx);

        TEMPLATE(T, poly_div_series) (c, a, b, n, ctx);
        TEMPLATE(T, poly_mullow) (d, c, b, n, ctx);

        result = (TEMPLATE(T, poly_equal_trunc) (a, d, n, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print) (a, ctx),
                flint_printf("\n\n");
            flint_printf("b = "), TEMPLATE(T, poly_print) (b, ctx),
                flint_printf("\n\n");
            flint_printf("c = "), TEMPLATE(T, poly_print) (c, ctx),
                flint_printf("\n\n");
            flint_printf("d = "), TEMPLATE(T, poly_print) (d, ctx),
                flint_printf("\n\n");
            flint_printf("ctx = "), TEMPLATE(T, ctx_print) (ctx),
                flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (d, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
