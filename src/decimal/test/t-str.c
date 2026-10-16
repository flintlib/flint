/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "decimal.h"
#include "fmpq.h"
#include "gr.h"

TEST_FUNCTION_START(decfloat_str, state)
{
    slong iter;

    /* fixed examples */
    {
        gr_ctx_t ctx;
        decfloat_t x;
        char * s;
        slong i;

        const char * inputs[] = { "0", "1", "-1", "1.5", "0.001", "123456789", "1e5", "1.5E-3",
            "-0.000000123", "1e21", "1e20", "123.456e-2", "0.1", "00012.500", "inf", "-inf", "nan",
            "1e-7", "1e-6", "12345678901234567890123", "-3.14159", NULL };
        const char * outputs[] = { "0", "1", "-1", "1.5", "0.001", "123456789", "100000", "0.0015",
            "-1.23e-7", "1e21", "100000000000000000000", "1.23456", "0.1", "12.5", "inf", "-inf", "nan",
            "1e-7", "0.000001", "1.2345678901234567890123e22", "-3.14159", NULL };

        gr_ctx_init_decfloat(ctx, DECIMAL_PREC_EXACT, DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN);
        decfloat_init(x, ctx);

        for (i = 0; inputs[i] != NULL; i++)
        {
            GR_MUST_SUCCEED(decfloat_set_str(x, inputs[i], ctx));
            s = decfloat_get_str(x, ctx);
            if (strcmp(s, outputs[i]))
            {
                flint_printf("FAIL: string example %wd\n", i);
                flint_printf("input = %s, output = %s, expected = %s\n", inputs[i], s, outputs[i]);
                flint_abort();
            }
            flint_free(s);
        }

        /* rounding when parsing */
        decimal_ctx_set_prec(ctx, 3);
        decimal_ctx_set_rnd(ctx, DECIMAL_RND_NEAR);
        GR_MUST_SUCCEED(decfloat_set_str(x, "12345", ctx));
        s = decfloat_get_str(x, ctx);
        if (strcmp(s, "12300")) { flint_printf("FAIL: rounding parse: %s\n", s); flint_abort(); }
        flint_free(s);
        GR_MUST_SUCCEED(decfloat_set_str(x, "0.12355", ctx));
        s = decfloat_get_str(x, ctx);
        if (strcmp(s, "0.124")) { flint_printf("FAIL: rounding parse: %s\n", s); flint_abort(); }
        flint_free(s);
        GR_MUST_SUCCEED(decfloat_set_str(x, "0.12345", ctx));
        s = decfloat_get_str(x, ctx);
        if (strcmp(s, "0.123")) { flint_printf("FAIL: rounding parse: %s\n", s); flint_abort(); }
        flint_free(s);
        GR_MUST_SUCCEED(decfloat_set_str(x, "2.5e-1000", ctx));
        s = decfloat_get_str(x, ctx);
        if (strcmp(s, "2.5e-1000")) { flint_printf("FAIL: parse: %s\n", s); flint_abort(); }
        flint_free(s);
        GR_MUST_SUCCEED(decfloat_set_str(x, "-1/3", ctx));
        s = decfloat_get_str(x, ctx);
        if (strcmp(s, "-0.333")) { flint_printf("FAIL: parse: %s\n", s); flint_abort(); }
        flint_free(s);

        /* invalid strings */
        if (decfloat_set_str(x, "1.2.3", ctx) == GR_SUCCESS) { flint_printf("FAIL: invalid string accepted\n"); flint_abort(); }
        if (decfloat_set_str(x, "abc", ctx) == GR_SUCCESS) { flint_printf("FAIL: invalid string accepted\n"); flint_abort(); }
        if (decfloat_set_str(x, "1e", ctx) == GR_SUCCESS) { flint_printf("FAIL: invalid string accepted\n"); flint_abort(); }

        decfloat_clear(x, ctx);
        gr_ctx_clear(ctx);
    }

    /* random roundtrips */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t x, y;
        char * s;
        int status;

        gr_ctx_init_decfloat_randtest(ctx, state, 60);
        decfloat_init(x, ctx);
        decfloat_init(y, ctx);

        GR_MUST_SUCCEED(decfloat_randtest(x, state, ctx));

        if (n_randint(state, 2))
            s = decfloat_get_str(x, ctx);
        else
            s = decfloat_get_str_sci(x, ctx);

        status = decfloat_set_round_str(y, s, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);

        if (status != GR_SUCCESS || (decfloat_equal(x, y, ctx) != T_TRUE && !DECFLOAT_IS_NAN(x)))
        {
            flint_printf("FAIL: string roundtrip\n");
            gr_ctx_println(ctx);
            flint_printf("x = %{gr}\n", x, ctx);
            flint_printf("s = %s\n", s);
            flint_printf("y = %{gr} (%d)\n", y, ctx, status);
            flint_abort();
        }

        flint_free(s);
        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
