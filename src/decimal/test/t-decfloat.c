/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "decimal.h"
#include "gr.h"

TEST_FUNCTION_START(decfloat, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;

        gr_ctx_init_decfloat_randtest(ctx, state, 60);

        if (gr_ctx_is_ring(ctx) == T_TRUE)
            gr_test_ring(ctx, 50, 0);
        else
            gr_test_floating_point(ctx, 50, 0);

        gr_ctx_clear(ctx);
    }

    /* the default context at a few standard precisions */
    {
        gr_ctx_t ctx;
        slong prec;

        for (prec = 1; prec <= 100; prec = 2 * prec + 1)
        {
            gr_ctx_init_decfloat(ctx, prec, DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN);
            gr_test_floating_point(ctx, 20 * flint_test_multiplier(), 0);
            gr_ctx_clear(ctx);

            gr_ctx_init_decfloat(ctx, prec, 0);
            gr_test_floating_point(ctx, 20 * flint_test_multiplier(), 0);
            gr_ctx_clear(ctx);
        }

        gr_ctx_init_decfloat(ctx, DECIMAL_PREC_EXACT, 0);
        gr_test_ring(ctx, 100 * flint_test_multiplier(), 0);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
