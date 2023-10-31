/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "qadic.h"

TEST_FUNCTION_START(qadic_trace, state)
{
    int i, result;

    /* Compare with sum of Galois conjugates */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c;
        padic_t x, y;
        slong j;
        int ans;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        padic_init2(x, N);
        padic_init2(y, N);

        qadic_randtest(a, state, ctx);

        qadic_trace(x, a, ctx);

        qadic_zero(b);
        for (j = 0; j < d; j++)
        {
            qadic_frobenius(c, a, j, ctx);
            qadic_add(b, b, c, ctx);
        }
        ans = qadic_get_padic(y, b, ctx);

        result = (ans && padic_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("x = "), padic_print(x, &ctx->pctx), flint_printf("\n");
            flint_printf("y = "), padic_print(y, &ctx->pctx), flint_printf("\n");
            for (j = 0; j < d; j++)
            {
                qadic_frobenius(c, a, j, ctx);
                flint_printf("sigma^%wd = ", j), qadic_print_pretty(c, ctx), flint_printf("\n");
            }
            flint_printf("ans = %d\n", ans);
            qadic_ctx_print(ctx);
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        padic_clear(x);
        padic_clear(y);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
