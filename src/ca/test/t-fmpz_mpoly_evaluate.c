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
#include "ca_vec.h"

TEST_FUNCTION_START(ca_fmpz_mpoly_evaluate, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t cactx;
        slong i, n;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        ca_ptr x;
        ca_t fx, gx, hx, y;

        ca_ctx_init(cactx);

        n = 1 + n_randint(state, 5);
        fmpz_mpoly_ctx_init(ctx, n, ORD_LEX);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        x = _ca_vec_init(n, cactx);
        ca_init(fx, cactx);
        ca_init(gx, cactx);
        ca_init(hx, cactx);
        ca_init(y, cactx);

        for (i = 0; i < n; i++)
            ca_randtest(x + i, state, 5, 5, cactx);

        fmpz_mpoly_randtest_bound(f, state, 1 + n_randint(state, 30), 10, 1 + n_randint(state, 4), ctx);
        fmpz_mpoly_randtest_bound(g, state, 1 + n_randint(state, 30), 10, 1 + n_randint(state, 4), ctx);
        fmpz_mpoly_add(h, f, g, ctx);

        ca_fmpz_mpoly_evaluate(fx, f, x, ctx, cactx);
        ca_fmpz_mpoly_evaluate(gx, g, x, ctx, cactx);
        ca_fmpz_mpoly_evaluate(hx, h, x, ctx, cactx);
        ca_add(y, fx, gx, cactx);

        if (ca_check_equal(y, hx, cactx) == T_FALSE)
        {
            flint_printf("FAIL!\n");
            flint_printf("f = "); fmpz_mpoly_print_pretty(f, NULL, ctx); flint_printf("\n\n");
            flint_printf("g = "); fmpz_mpoly_print_pretty(g, NULL, ctx); flint_printf("\n\n");
            flint_printf("h = "); fmpz_mpoly_print_pretty(h, NULL, ctx); flint_printf("\n\n");
            for (i = 0; i < n; i++)
            {
                flint_printf("x%wd = ", i + 1); ca_print(x + i, cactx); flint_printf("\n\n");
            }
            flint_printf("fx = "); ca_print(fx, cactx); flint_printf("\n\n");
            flint_printf("gx = "); ca_print(gx, cactx); flint_printf("\n\n");
            flint_printf("hx = "); ca_print(hx, cactx); flint_printf("\n\n");
            flint_printf("y = "); ca_print(y, cactx); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);

        _ca_vec_clear(x, n, cactx);
        ca_clear(fx, cactx);
        ca_clear(gx, cactx);
        ca_clear(hx, cactx);
        ca_clear(y, cactx);

        ca_ctx_clear(cactx);
    }

    TEST_FUNCTION_END(state);
}
