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

TEST_FUNCTION_START(ca_pow, state)
{
    slong iter;

    /* check numerical evaluation */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y, z;
        acb_t ax, ay, az, axy;
        slong prec;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(z, ctx);

        acb_init(ax);
        acb_init(ay);
        acb_init(az);
        acb_init(axy);

        prec = 10 + n_randint(state, 100);
        ca_randtest(x, state, 5, 5, ctx);
        ca_randtest(y, state, 5, 5, ctx);
        if (n_randint(state, 2))
            ca_exp(x, x, ctx);
        if (n_randint(state, 2))
        {
            ca_pow(x, x, y, ctx);
            ca_randtest(y, state, 5, 5, ctx);
        }

        ca_pow(z, x, y, ctx);

        ca_get_acb(ax, x, prec, ctx);
        ca_get_acb(ay, y, prec, ctx);
        ca_get_acb(az, z, prec, ctx);

        acb_pow(axy, ax, ay, prec);

        if (!acb_overlaps(axy, az))
        {
            flint_printf("FAIL (overlap)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); flint_printf("\n\n");
            flint_printf("z = "); ca_print(z, ctx); flint_printf("\n\n");
            flint_printf("ax = "); acb_printn(ax, 30, ARB_STR_NO_RADIUS); flint_printf("\n\n");
            flint_printf("ay = "); acb_printn(ay, 30, ARB_STR_NO_RADIUS); flint_printf("\n\n");
            flint_printf("az = "); acb_printn(az, 30, ARB_STR_NO_RADIUS); flint_printf("\n\n");
            flint_printf("axy = "); acb_printn(axy, 30, ARB_STR_NO_RADIUS); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(ax);
        acb_clear(ay);
        acb_clear(az);
        acb_clear(axy);

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(z, ctx);
        ca_ctx_clear(ctx);
    }

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, a, b, xa, xb, xaxb, ab, xab;
        truth_t equal;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(a, ctx);
        ca_init(b, ctx);
        ca_init(xa, ctx);
        ca_init(xb, ctx);
        ca_init(xaxb, ctx);
        ca_init(ab, ctx);
        ca_init(xab, ctx);

        /* otherwise this will be too slow */
        ctx->options[CA_OPT_QQBAR_DEG_LIMIT] = 40;

        /* x^a * x^b = x^(a+b) */

        do {
            ca_randtest(x, state, 5, 5, ctx);
        } while (ca_check_is_zero(x, ctx) != T_FALSE);

        ca_randtest(a, state, 5, 5, ctx);
        ca_randtest(b, state, 5, 5, ctx);

        ca_pow(xa, x, a, ctx);
        ca_pow(xb, x, b, ctx);
        ca_mul(xaxb, xa, xb, ctx);

        ca_add(ab, a, b, ctx);
        ca_pow(xab, x, ab, ctx);

        equal = ca_check_equal(xab, xaxb, ctx);

        if (equal == T_FALSE)
        {
            flint_printf("FAIL x^a * x^b = x^(a+b)\n\n");
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n");
            flint_printf("a = "); ca_print(a, ctx); flint_printf("\n\n");
            flint_printf("b = "); ca_print(b, ctx); flint_printf("\n\n");
            flint_printf("xa = "); ca_print(xa, ctx); flint_printf("\n\n");
            flint_printf("xb = "); ca_print(xb, ctx); flint_printf("\n\n");
            flint_printf("xaxb = "); ca_print(xaxb, ctx); flint_printf("\n\n");
            flint_printf("ab = "); ca_print(ab, ctx); flint_printf("\n\n");
            flint_printf("xab = "); ca_print(xab, ctx); flint_printf("\n\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(a, ctx);
        ca_clear(b, ctx);
        ca_clear(xa, ctx);
        ca_clear(xb, ctx);
        ca_clear(xaxb, ctx);
        ca_clear(ab, ctx);
        ca_clear(xab, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
