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

TEST_FUNCTION_START(ca_log, state)
{
    slong iter;

    /* check numerical evaluation */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y;
        acb_t ax, ay, logax;
        slong prec;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);

        acb_init(ax);
        acb_init(ay);
        acb_init(logax);

        prec = 10 + n_randint(state, 100);
        ca_randtest(x, state, 5, 5, ctx);
        ca_randtest(y, state, 5, 5, ctx);
        if (n_randint(state, 2))
            ca_exp(x, x, ctx);
        if (n_randint(state, 2))
            ca_pow(x, x, y, ctx);

        ca_log(y, x, ctx);

        ca_get_acb(ax, x, prec, ctx);
        ca_get_acb(ay, y, prec, ctx);
        acb_log(logax, ax, prec);

        if (!acb_overlaps(logax, ay))
        {
            flint_printf("FAIL (overlap)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_printf("ax = "); acb_printn(ax, 30, ARB_STR_NO_RADIUS); flint_printf("\n\n");
            flint_printf("ay = "); acb_printn(ay, 30, ARB_STR_NO_RADIUS); flint_printf("\n\n");
            flint_printf("logax = "); acb_printn(logax, 30, ARB_STR_NO_RADIUS); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(ax);
        acb_clear(ay);
        acb_clear(logax);

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_ctx_clear(ctx);
    }

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y, z, a, b, c, d, e, f;
        truth_t equal, zero;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(z, ctx);
        ca_init(a, ctx);
        ca_init(b, ctx);
        ca_init(c, ctx);
        ca_init(d, ctx);
        ca_init(e, ctx);
        ca_init(f, ctx);

        /* log(x * y * z) - log(x) - log(y) - log(z) */
        /* testing for positive x and y here -- need better test code for complex numbers */

        ca_randtest(x, state, 5, 5, ctx);
        ca_randtest(y, state, 5, 5, ctx);
        ca_randtest(z, state, 5, 5, ctx);

        ca_abs(x, x, ctx);
        ca_abs(y, y, ctx);

        /* Test a bug */
        if (iter == 0)
        {
            ca_sqrt_ui(x, 2, ctx);
            ca_div_ui(x, x, 3, ctx);
            ca_sqrt(x, x, ctx);

            ca_set_ui(y, 1, ctx);
            ca_div_ui(y, y, 2, ctx);
            ca_one(z, ctx);
        }

        /* a = log(x*y*z) */
        ca_mul(a, x, y, ctx);
        ca_mul(a, a, z, ctx);
        ca_log(a, a, ctx);

        /* b = log(x), c = log(y), d = log(z) */
        ca_log(b, x, ctx);
        ca_log(c, y, ctx);
        ca_log(d, z, ctx);

        /* e = log(x) + log(y) + log(c) */
        ca_add(e, b, c, ctx);
        ca_add(e, e, d, ctx);

        ca_sub(f, a, e, ctx);

        equal = ca_check_equal(a, e, ctx);

        if (ca_check_is_infinity(a, ctx) == T_FALSE)
        {
            zero = ca_check_is_zero(f, ctx);
        }
        else
        {
            zero = T_UNKNOWN;
        }

        if (equal == T_FALSE || zero == T_FALSE)
        {
            flint_printf("FAIL (log(x * y * z) - log(x) - log(y) - log(z) != 0)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_printf("z = "); ca_print(z, ctx); flint_printf("\n\n");
            flint_printf("a = "); ca_print(a, ctx); flint_printf("\n\n");
            flint_printf("b = "); ca_print(b, ctx); flint_printf("\n\n");
            flint_printf("c = "); ca_print(c, ctx); flint_printf("\n\n");
            flint_printf("d = "); ca_print(d, ctx); flint_printf("\n\n");
            flint_printf("e = "); ca_print(e, ctx); flint_printf("\n\n");
            flint_printf("f = "); ca_print(f, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(z, ctx);
        ca_clear(a, ctx);
        ca_clear(b, ctx);
        ca_clear(c, ctx);
        ca_clear(d, ctx);
        ca_clear(e, ctx);
        ca_clear(f, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
