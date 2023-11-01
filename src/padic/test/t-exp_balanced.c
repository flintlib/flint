/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "padic.h"

TEST_FUNCTION_START(padic_exp_balanced, state)
{
    int i, result;

/** p == 2 *******************************************************************/

    /* Check aliasing: a = exp(a) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b;
        int ans1, ans2;

        fmpz_init_set_ui(p, 2);

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;

        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);

        padic_randtest(a, state, ctx);

        ans1 = padic_exp_balanced(b, a, ctx);
        ans2 = padic_exp_balanced(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, b)));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("ans1 = %d\n", ans1);
            flint_printf("ans2 = %d\n", ans2);
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Functional equation: exp(a + b) == exp(a) exp(b) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;
        int ans1, ans2, ans3;

        fmpz_init_set_ui(p, 2);

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;

        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);
        padic_init2(e, N);
        padic_init2(f, N);
        padic_init2(g, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);
        padic_add(c, a, b, ctx);

        ans1 = padic_exp_balanced(d, a, ctx);
        ans2 = padic_exp_balanced(e, b, ctx);
        padic_mul(f, d, e, ctx);

        ans3 = padic_exp_balanced(g, c, ctx);

        result = (!ans1 || !ans2 || (ans3 && padic_equal(f, g)));
        if (!result)
        {
            flint_printf("FAIL (functional equation):\n\n");
            flint_printf("a                 = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b                 = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("c = a + b         = "), padic_print(c, ctx), flint_printf("\n");
            flint_printf("d = exp(a)        = "), padic_print(d, ctx), flint_printf("\n");
            flint_printf("e = exp(b)        = "), padic_print(e, ctx), flint_printf("\n");
            flint_printf("f = exp(a) exp(b) = "), padic_print(f, ctx), flint_printf("\n");
            flint_printf("g = exp(a + b)    = "), padic_print(g, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);
        padic_clear(e);
        padic_clear(f);
        padic_clear(g);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

/** p > 2 ********************************************************************/

    /* Check aliasing: a = exp(a) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b;
        int ans1, ans2;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;

        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);

        padic_randtest(a, state, ctx);

        ans1 = padic_exp_balanced(b, a, ctx);
        ans2 = padic_exp_balanced(a, a, ctx);

        result = ((ans1 == ans2) && (!ans1 || padic_equal(a, b)));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("ans1 = %d\n", ans1);
            flint_printf("ans2 = %d\n", ans2);
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Functional equation: exp(a + b) == exp(a) exp(b) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e, f, g;
        int ans1, ans2, ans3;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;

        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);
        padic_init2(e, N);
        padic_init2(f, N);
        padic_init2(g, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);
        padic_add(c, a, b, ctx);

        ans1 = padic_exp_balanced(d, a, ctx);
        ans2 = padic_exp_balanced(e, b, ctx);
        padic_mul(f, d, e, ctx);

        ans3 = padic_exp_balanced(g, c, ctx);

        result = (!ans1 || !ans2 || (ans3 && padic_equal(f, g)));
        if (!result)
        {
            flint_printf("FAIL (functional equation):\n\n");
            flint_printf("a                 = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b                 = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("c = a + b         = "), padic_print(c, ctx), flint_printf("\n");
            flint_printf("d = exp(a)        = "), padic_print(d, ctx), flint_printf("\n");
            flint_printf("e = exp(b)        = "), padic_print(e, ctx), flint_printf("\n");
            flint_printf("f = exp(a) exp(b) = "), padic_print(f, ctx), flint_printf("\n");
            flint_printf("g = exp(a + b)    = "), padic_print(g, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);
        padic_clear(e);
        padic_clear(f);
        padic_clear(g);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
