/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca.h"

TEST_FUNCTION_START(ca_re_im, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, rex, imx, y, absx;
        truth_t equal;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(rex, ctx);
        ca_init(imx, ctx);
        ca_init(y, ctx);
        ca_init(absx, ctx);

        ca_randtest(x, state, 5, 5, ctx);

        ca_re(rex, x, ctx);
        ca_im(imx, x, ctx);

        /* test re(x) + im(x)*i = x */
        ca_i(y, ctx);
        ca_mul(y, y, imx, ctx);
        ca_add(y, y, rex, ctx);

        equal = ca_check_equal(x, y, ctx);

        if (equal == T_FALSE)
        {
            flint_printf("FAIL (re(x) + im(x)*i != x)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("rex = "); ca_print(rex, ctx); flint_printf("\n\n");
            flint_printf("imx = "); ca_print(imx, ctx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_abort();
        }

        CA_TEST_PROPERTY(ca_check_is_real, "is_real", rex, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", imx, ctx, T_TRUE);

        /* test sqrt(re(x)^2 + im(x)^2) = abs(x) */
        ca_sqr(y, rex, ctx);
        ca_sqr(absx, imx, ctx);
        ca_add(y, y, absx, ctx);
        ca_sqrt(y, y, ctx);
        ca_abs(absx, x, ctx);

        equal = ca_check_equal(absx, y, ctx);

        if (equal == T_FALSE)
        {
            flint_printf("FAIL (sqrt(re(x)^2 + im(x)^2) != abs(x))\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("rex = "); ca_print(rex, ctx); flint_printf("\n\n");
            flint_printf("imx = "); ca_print(imx, ctx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_printf("absx = "); ca_print(absx, ctx); flint_printf("\n\n");
            flint_abort();
        }

        CA_TEST_PROPERTY(ca_check_is_real, "is_real", y, ctx, T_TRUE);
        CA_TEST_PROPERTY(ca_check_is_real, "is_real", absx, ctx, T_TRUE);

        ca_clear(x, ctx);
        ca_clear(rex, ctx);
        ca_clear(imx, ctx);
        ca_clear(y, ctx);
        ca_clear(absx, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
