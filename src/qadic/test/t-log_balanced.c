/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "qadic.h"

TEST_FUNCTION_START(qadic_log_balanced, state)
{
    int i, result;

/** p == 2 *******************************************************************/

/** p > 2 ********************************************************************/

    /* Check aliasing: a = log(a) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        slong ans1, ans2;
        qadic_ctx_t ctx;

        qadic_t a, b;

        fmpz_init_set_ui(p, n_randprime(state, 3 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);

        qadic_randtest_not_zero(a, state, ctx);
        if (a->val < 1)
            a->val = 1;
        padic_poly_reduce(a, &ctx->pctx);
        qadic_one(b);
        qadic_add(a, a, b, ctx);

        ans1 = qadic_log_balanced(b, a, ctx);
        ans2 = qadic_log_balanced(a, a, ctx);

        result = (ans1 == ans2) && (!ans1 || qadic_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
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

    /* Check: log(a) + log(b) == log(a * b) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong deg, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c, d, e, f, g;

        fmpz_init_set_ui(p, n_randprime(state, 3 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        qadic_init2(d, N);
        qadic_init2(e, N);
        qadic_init2(f, N);
        qadic_init2(g, N);

        qadic_randtest_not_zero(a, state, ctx);
        if (a->val < 1)
            a->val = 1;
        padic_poly_reduce(a, &ctx->pctx);
        qadic_one(c);
        qadic_add(a, a, c, ctx);

        qadic_randtest_not_zero(b, state, ctx);
        if (b->val < 1)
            b->val = 1;
        padic_poly_reduce(b, &ctx->pctx);
        qadic_one(c);
        qadic_add(b, b, c, ctx);

        qadic_mul(c, a, b, ctx);

        qadic_log_balanced(d, a, ctx);
        qadic_log_balanced(e, b, ctx);
        qadic_add(f, d, e, ctx);

        qadic_log_balanced(g, c, ctx);

        result = (qadic_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL (functional equation):\n\n");
            flint_printf("a                   = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b                   = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = a * b           = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("d = log(a)          = "), qadic_print_pretty(d, ctx), flint_printf("\n");
            flint_printf("e = log(b)          = "), qadic_print_pretty(e, ctx), flint_printf("\n");
            flint_printf("f = log(a) + log(b) = "), qadic_print_pretty(f, ctx), flint_printf("\n");
            flint_printf("g = log(a * b)      = "), qadic_print_pretty(g, ctx), flint_printf("\n");
            qadic_ctx_print(ctx);
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);
        qadic_clear(d);
        qadic_clear(e);
        qadic_clear(f);
        qadic_clear(g);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Check: log(exp(x)) == x */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong deg, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c;

        fmpz_init_set_ui(p, n_randprime(state, 3 + n_randint(state, 3), 1));
        deg = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, deg, FLINT_MAX(0,N-10), FLINT_MAX(0,N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);

        qadic_randtest_val(a, state, (*p == WORD(2)) + 1, ctx);

        qadic_exp(b, a, ctx);
        qadic_log_balanced(c, b, ctx);

        result = (qadic_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (log(exp(x)) == x):\n\n");
            flint_printf("a =          "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = exp(a) = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = log(b) = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
