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

TEST_FUNCTION_START(ca_neg, state)
{
    {
        ca_ctx_t ctx;
        ca_t x, y, z;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        ca_init(y, ctx);
        ca_init(z, ctx);

        ca_pos_inf(x, ctx);
        ca_neg(y, x, ctx);
        ca_set(z, x, ctx);
        ca_neg(z, z, ctx);

        if (ca_check_is_neg_inf(y, ctx) != T_TRUE || ca_check_is_neg_inf(z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: pos_inf\n");
            flint_abort();
        }

        ca_neg_inf(x, ctx);
        ca_neg(y, x, ctx);
        ca_set(z, x, ctx);
        ca_neg(z, z, ctx);

        if (ca_check_is_pos_inf(y, ctx) != T_TRUE || ca_check_is_pos_inf(z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: neg_inf\n");
            flint_abort();
        }

        ca_pos_i_inf(x, ctx);
        ca_neg(y, x, ctx);
        ca_set(z, x, ctx);
        ca_neg(z, z, ctx);

        if (ca_check_is_neg_i_inf(y, ctx) != T_TRUE || ca_check_is_neg_i_inf(z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: pos_i_inf\n");
            flint_abort();
        }

        ca_undefined(x, ctx);
        ca_neg(y, x, ctx);
        ca_set(z, x, ctx);
        ca_neg(z, z, ctx);

        if (ca_check_is_undefined(y, ctx) != T_TRUE || ca_check_is_undefined(z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: neg_i_inf\n");
            flint_abort();
        }

        ca_neg_i_inf(x, ctx);
        ca_neg(y, x, ctx);
        ca_set(z, x, ctx);
        ca_neg(z, z, ctx);

        if (ca_check_is_pos_i_inf(y, ctx) != T_TRUE || ca_check_is_pos_i_inf(z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: neg_i_inf\n");
            flint_abort();
        }

        ca_uinf(x, ctx);
        ca_neg(y, x, ctx);
        ca_set(z, x, ctx);
        ca_neg(z, z, ctx);

        if (ca_check_is_uinf(y, ctx) != T_TRUE || ca_check_is_uinf(z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: neg_i_inf\n");
            flint_abort();
        }

        ca_i(x, ctx);
        ca_neg(y, x, ctx);
        ca_set(z, x, ctx);
        ca_neg(z, z, ctx);

        if (ca_check_is_neg_i(y, ctx) != T_TRUE || ca_check_is_neg_i(z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: i\n");
            flint_abort();
        }

        ca_one(x, ctx);
        ca_neg(y, x, ctx);
        ca_set(z, x, ctx);
        ca_neg(z, z, ctx);

        if (ca_check_is_neg_one(y, ctx) != T_TRUE || ca_check_is_neg_one(z, ctx) != T_TRUE)
        {
            flint_printf("FAIL: one\n");
            flint_abort();
        }

        ca_unknown(x, ctx);
        ca_neg(y, x, ctx);
        ca_set(z, x, ctx);
        ca_neg(z, z, ctx);

        if (!ca_is_unknown(y, ctx) || !ca_is_unknown(z, ctx))
        {
            flint_printf("FAIL: unknown\n");
            flint_abort();
        }

        ca_clear(x, ctx);
        ca_clear(y, ctx);
        ca_clear(z, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
