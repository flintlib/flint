/*
    Copyright (C) 2009 William Hart
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

TEST_TEMPLATE_FUNCTION_START(T, poly_inv_series_newton, state)
{
    int i, result;

    /* Check Q^{-1} * Q is congruent 1 mod t^n */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) a, b, c, one;
        slong n = n_randint(state, 80) + 1;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (c, ctx);
        TEMPLATE(T, poly_init) (one, ctx);

        TEMPLATE(T, poly_randtest_not_zero) (a, state,
                                             n_randint(state, 80) + 1, ctx);
        TEMPLATE(T, randtest_not_zero) (a->coeffs, state, ctx);
        TEMPLATE(T, poly_one) (one, ctx);

        TEMPLATE(T, poly_inv_series_newton) (b, a, n, ctx);
        TEMPLATE(T, poly_mullow) (c, a, b, n, ctx);

        result = (TEMPLATE(T, poly_equal) (c, one, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), TEMPLATE(T, poly_print) (a, ctx),
                flint_printf("\n\n");
            flint_printf("b = "), TEMPLATE(T, poly_print) (b, ctx),
                flint_printf("\n\n");
            flint_printf("c = "), TEMPLATE(T, poly_print) (c, ctx),
                flint_printf("\n\n");
            flint_printf("ctx = "), TEMPLATE(T, ctx_print) (ctx),
                flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (c, ctx);
        TEMPLATE(T, poly_clear) (one, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
