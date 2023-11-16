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

void
ca_randtest_zero(ca_t x, flint_rand_t state, ca_ctx_t ctx)
{
    if (n_randint(state, 2))
    {
        ca_zero(x, ctx);
    }
    else
    {
        ca_t t;
        ca_init(t, ctx);
        ca_sqrt_ui(x, 6, ctx);
        ca_mul_ui(x, x, 2, ctx);
        ca_add_ui(x, x, 5, ctx);
        ca_sqrt(x, x, ctx);
        ca_sqrt_ui(t, 2, ctx);
        ca_sub(x, x, t, ctx);
        ca_sqrt_ui(t, 3, ctx);
        ca_sub(x, x, t, ctx);
        ca_clear(t, ctx);
    }
}

TEST_FUNCTION_START(ca_erf, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x, y, v, s;
        acb_t ax, ay, as, as2;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(v, ctx);
        ca_init(s, ctx);
        acb_init(ax);
        acb_init(ay);
        acb_init(as);
        acb_init(as2);

        ca_randtest(x, state, 3, 5, ctx);
        ca_randtest_zero(v, state, ctx);
        ca_add(y, x, v, ctx);

        switch (n_randint(state, 4))
        {
            case 0: ca_neg(y, y, ctx); break;
            case 1: ca_i(s, ctx); ca_mul(y, y, s, ctx); break;
            case 2: ca_neg_i(s, ctx); ca_mul(y, y, s, ctx); break;
        }

        switch (n_randint(state, 3))
        {
            case 0: ca_erf(x, x, ctx); break;
            case 1: ca_erfc(x, x, ctx); break;
            default: ca_erfi(x, x, ctx); break;
        }

        switch (n_randint(state, 3))
        {
            case 0: ca_erf(y, y, ctx); break;
            case 1: ca_erfc(y, y, ctx); break;
            default: ca_erfi(y, y, ctx); break;
        }

        if (n_randint(state, 2))
            ca_sqr(x, x, ctx);

        switch (n_randint(state, 4))
        {
            case 0: ca_neg(x, x, ctx); break;
            case 1: ca_i(s, ctx); ca_mul(x, x, s, ctx); break;
            case 2: ca_neg_i(s, ctx); ca_mul(x, x, s, ctx); break;
        }

        if (n_randint(state, 2))
            ca_sqr(y, y, ctx);

        switch (n_randint(state, 4))
        {
            case 0: ca_neg(y, y, ctx); break;
            case 1: ca_i(s, ctx); ca_mul(y, y, s, ctx); break;
            case 2: ca_neg_i(s, ctx); ca_mul(y, y, s, ctx); break;
        }

        ca_add(s, x, y, ctx);

        ca_get_acb(ax, x, 53, ctx);
        ca_get_acb(ay, y, 53, ctx);
        ca_get_acb(as, s, 53, ctx);

        acb_add(as2, ax, ay, 53);

        if (!acb_overlaps(as, as2))
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); ca_print(x, ctx); printf("\n\n");
            flint_printf("y = "); ca_print(y, ctx); printf("\n\n");
            flint_printf("s = "); ca_print(s, ctx); printf("\n\n");
            flint_abort();
        }

        acb_clear(ax);
        acb_clear(ay);
        acb_clear(as);
        acb_clear(as2);
        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(v, ctx);
        ca_clear(s, ctx);
        ca_ctx_clear(ctx);
    }

    {
        ca_ctx_t ctx;
        ca_t x, y;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);

        /* Erf(2*Log(Sqrt(1/2-Sqrt(2)/4))+Log(4)) - Erf(Log(2-Sqrt(2))) */
        ca_sqrt_ui(x, 2, ctx);
        ca_div_ui(x, x, 4, ctx);
        ca_one(y, ctx);
        ca_div_ui(y, y, 2, ctx);
        ca_sub(x, y, x, ctx);
        ca_sqrt(x, x, ctx);
        ca_log(x, x, ctx);
        ca_mul_ui(x, x, 2, ctx);
        ca_set_ui(y, 4, ctx);
        ca_log(y, y, ctx);
        ca_add(x, y, x, ctx);
        ca_erf(x, x, ctx);
        ca_sqrt_ui(y, 2, ctx);
        ca_ui_sub(y, 2, y, ctx);
        ca_log(y, y, ctx);
        ca_erf(y, y, ctx);
        ca_sub(x, x, y, ctx);

        if (ca_check_is_zero(x, ctx) != T_TRUE)
        {
            flint_printf("FAIL (example 1)\n");
            flint_printf("x = "); ca_print(x, ctx); printf("\n\n");
            flint_abort();
        }

        /* Erf(Sqrt(2))^2 + Erfi(Sqrt(-2))^2 */
        ca_sqrt_ui(x, 2, ctx);
        ca_erf(x, x, ctx);
        ca_pow_ui(x, x, 2, ctx);
        ca_set_si(y, -2, ctx);
        ca_sqrt(y, y, ctx);
        ca_erfi(y, y, ctx);
        ca_pow_ui(y, y, 2, ctx);
        ca_add(x, x, y, ctx);

        if (ca_check_is_zero(x, ctx) != T_TRUE)
        {
            flint_printf("FAIL (example 2)\n");
            flint_printf("x = "); ca_print(x, ctx); printf("\n\n");
            flint_abort();
        }

        /* Erf(Sqrt(2))^4 - Erfi(Sqrt(-2))^4 */
        ca_sqrt_ui(x, 2, ctx);
        ca_erf(x, x, ctx);
        ca_pow_ui(x, x, 4, ctx);
        ca_set_si(y, -2, ctx);
        ca_sqrt(y, y, ctx);
        ca_erfi(y, y, ctx);
        ca_pow_ui(y, y, 4, ctx);
        ca_sub(x, x, y, ctx);

        if (ca_check_is_zero(x, ctx) != T_TRUE)
        {
            flint_printf("FAIL (example 3)\n");
            flint_printf("x = "); ca_print(x, ctx); printf("\n\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
