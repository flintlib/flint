/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_evaluate_fmpz_mpoly, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong i, n;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        qqbar_ptr x;
        qqbar_t fx, gx, hx, y;
        int s1, s2, s3;

        n = 1 + n_randint(state, 5);
        fmpz_mpoly_ctx_init(ctx, n, ORD_LEX);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        x = _qqbar_vec_init(n);
        qqbar_init(fx);
        qqbar_init(gx);
        qqbar_init(hx);
        qqbar_init(y);

        for (i = 0; i < n; i++)
            qqbar_randtest(x + i, state, 1 + n_randint(state, 2), 10);

        fmpz_mpoly_randtest_bound(f, state, 1 + n_randint(state, 30), 10, 1 + n_randint(state, 4), ctx);
        fmpz_mpoly_randtest_bound(g, state, 1 + n_randint(state, 30), 10, 1 + n_randint(state, 4), ctx);
        fmpz_mpoly_add(h, f, g, ctx);

        s1 = qqbar_evaluate_fmpz_mpoly(fx, f, x, 40, 1000, ctx);
        s2 = qqbar_evaluate_fmpz_mpoly(gx, g, x, 40, 1000, ctx);
        s3 = qqbar_evaluate_fmpz_mpoly(hx, h, x, 40, 1000, ctx);
        qqbar_add(y, fx, gx);

        if (s1 && s2 && s3 && !qqbar_equal(y, hx))
        {
            flint_printf("FAIL!\n");
            flint_printf("f = "); fmpz_mpoly_print_pretty(f, NULL, ctx); flint_printf("\n\n");
            flint_printf("g = "); fmpz_mpoly_print_pretty(g, NULL, ctx); flint_printf("\n\n");
            flint_printf("h = "); fmpz_mpoly_print_pretty(h, NULL, ctx); flint_printf("\n\n");
            for (i = 0; i < n; i++)
            {
                flint_printf("x%wd = ", i + 1); qqbar_print(x + i); flint_printf("\n\n");
            }
            flint_printf("fx = "); qqbar_print(fx); flint_printf("\n\n");
            flint_printf("gx = "); qqbar_print(gx); flint_printf("\n\n");
            flint_printf("hx = "); qqbar_print(hx); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);

        _qqbar_vec_clear(x, n);
        qqbar_clear(fx);
        qqbar_clear(gx);
        qqbar_clear(hx);
        qqbar_clear(y);
    }

    TEST_FUNCTION_END(state);
}
