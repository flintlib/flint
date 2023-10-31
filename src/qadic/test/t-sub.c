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

TEST_FUNCTION_START(qadic_sub, state)
{
    int i, result;

    /* Check aliasing: a = a - b */
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

        qadic_randtest(a, state, ctx);
        qadic_randtest(b, state, ctx);

        qadic_sub(c, a, b, ctx);
        qadic_sub(a, a, b, ctx);

        result = (qadic_equal(a, c));
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

    /* Check aliasing: b = a - b */
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

        qadic_randtest(a, state, ctx);
        qadic_randtest(b, state, ctx);

        qadic_sub(c, a, b, ctx);
        qadic_sub(b, a, b, ctx);

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

    /* Check aliasing: a = a - a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, c;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(c, N);

        qadic_randtest(a, state, ctx);

        qadic_sub(c, a, a, ctx);
        qadic_sub(a, a, a, ctx);

        result = (qadic_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Check that a - b == -(b - a) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c1, c2;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c1, N);
        qadic_init2(c2, N);

        qadic_randtest(a, state, ctx);
        qadic_randtest(b, state, ctx);

        qadic_sub(c1, a, b, ctx);
        qadic_sub(c2, b, a, ctx);
        qadic_neg(c2, c2, ctx);

        result = (qadic_equal(c1, c2));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a  = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b  = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c1 = "), qadic_print_pretty(c1, ctx), flint_printf("\n");
            flint_printf("c2 = "), qadic_print_pretty(c2, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c1);
        qadic_clear(c2);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Check that (a - b) - c == a - (b + c) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c, lhs, rhs;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        qadic_init2(lhs, N);
        qadic_init2(rhs, N);

        qadic_randtest(a, state, ctx);
        qadic_randtest(b, state, ctx);
        qadic_randtest(c, state, ctx);

        qadic_sub(lhs, a, b, ctx);
        qadic_sub(lhs, lhs, c, ctx);
        qadic_add(rhs, b, c, ctx);
        qadic_sub(rhs, a, rhs, ctx);

        result = (qadic_equal(lhs, rhs));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a   = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b   = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c   = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("lhs = "), qadic_print_pretty(lhs, ctx), flint_printf("\n");
            flint_printf("rhs = "), qadic_print_pretty(rhs, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        qadic_clear(lhs);
        qadic_clear(rhs);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
