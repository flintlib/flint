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

TEST_FUNCTION_START(qadic_inv, state)
{
    int i, result;

    /* Check aliasing: a = ~a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);

        qadic_randtest_not_zero(a, state, ctx);
        qadic_set(b, a, ctx);

        qadic_inv(c, b, ctx);
        qadic_inv(b, b, ctx);

        result = (qadic_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Check a * ~a == 1 for units */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);

        do qadic_randtest_val(a, state, 0, ctx);
        while (qadic_is_zero(a));

        qadic_inv(b, a, ctx);
        qadic_mul(b, a, b, ctx);

        result = qadic_is_one(b);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
