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

TEST_FUNCTION_START(qadic_pow, state)
{
    int i, result;

    /* Check aliasing: a = a^e */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b;
        fmpz_t e;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        fmpz_init(e);

        qadic_randtest(a, state, ctx);
        fmpz_randtest_unsigned(e, state, 6);

        qadic_pow(b, a, e, ctx);
        qadic_pow(a, a, e, ctx);

        result = (qadic_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL (alias):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        fmpz_clear(e);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Compare with multiplication, for integral values */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c;
        fmpz_t e, f;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        fmpz_init(f);
        fmpz_init(e);

        qadic_randtest_int(a, state, ctx);
        fmpz_randtest_unsigned(e, state, 6);

        qadic_pow(b, a, e, ctx);
        qadic_one(c);
        for (fmpz_one(f); fmpz_cmp(f, e) <= 0; fmpz_add_ui(f, f, 1))
        {
            qadic_mul(c, c, a, ctx);
        }

        result = (qadic_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (cmp with mul):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("e = "), fmpz_print(e), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        fmpz_clear(e);
        fmpz_clear(f);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
