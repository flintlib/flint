/*
    Copyright (C) 2012, 2013 Sebastian Pancratz

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

TEST_FUNCTION_START(qadic_frobenius, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c;
        slong e;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = z_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);

        qadic_randtest(a, state, ctx);
        qadic_set(b, a, ctx);
        e = n_randint(state, 10) % d;

        qadic_frobenius(c, b, e, ctx);
        qadic_frobenius(b, b, e, ctx);

        result = (qadic_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (alias):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("e = %wd\n", e);
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(c);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Check sigma^e(x) == x^{p^e} mod p for integral values */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, c, lhs, rhs;
        slong e;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 50) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(c, N);
        qadic_init2(lhs, 1);
        qadic_init2(rhs, 1);

        qadic_randtest_int(a, state, ctx);
        e = n_randint(state, 10) % d;

        qadic_frobenius(b, a, e, ctx);
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, p, e);
            qadic_pow(c, a, t, ctx);
            fmpz_clear(t);
        }

        qadic_set(lhs, b, ctx);
        qadic_set(rhs, c, ctx);

        result = (qadic_equal(lhs, rhs));
        if (!result)
        {
            flint_printf("FAIL (sigma^e(x) = x^{p^e} mod p):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), qadic_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("lhs = "), qadic_print_pretty(lhs, ctx), flint_printf("\n");
            flint_printf("rhs = "), qadic_print_pretty(rhs, ctx), flint_printf("\n");
            flint_printf("e = %wd\n", e);
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

    /* Check sigma^e(x + y) = sigma^e(x) + sigma^e(y) on Zq */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, s, s1, s2, lhs, rhs;
        slong e;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 10) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(s, N);
        qadic_init2(s1, N);
        qadic_init2(s2, N);
        qadic_init2(lhs, N);
        qadic_init2(rhs, N);

        qadic_randtest_int(a, state, ctx);
        qadic_randtest_int(b, state, ctx);
        e = n_randint(state, 10) % d;

        qadic_add(s, a, b, ctx);
        qadic_frobenius(lhs, s, e, ctx);
        qadic_frobenius(s1, a, e, ctx);
        qadic_frobenius(s2, b, e, ctx);
        qadic_add(rhs, s1, s2, ctx);

        result = (qadic_equal(lhs, rhs));
        if (!result)
        {
            flint_printf("FAIL (sigma(a+b) = sigma(a) + sigma(b)):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("s = "), qadic_print_pretty(s, ctx), flint_printf("\n");
            flint_printf("s1 = "), qadic_print_pretty(s1, ctx), flint_printf("\n");
            flint_printf("s2 = "), qadic_print_pretty(s2, ctx), flint_printf("\n");
            flint_printf("lhs = "), qadic_print_pretty(lhs, ctx), flint_printf("\n");
            flint_printf("rhs = "), qadic_print_pretty(rhs, ctx), flint_printf("\n");
            flint_printf("e = %wd\n", e);
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(s);
        qadic_clear(s1);
        qadic_clear(s2);
        qadic_clear(lhs);
        qadic_clear(rhs);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    /* Check sigma^e(x * y) = sigma^e(x) * sigma^e(y) on Zq */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d, N;
        qadic_ctx_t ctx;

        qadic_t a, b, s, s1, s2, lhs, rhs;
        slong e;

        fmpz_init_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        N = n_randint(state, 10) + 1;
        qadic_ctx_init_conway(ctx, p, d, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), "a", PADIC_SERIES);

        qadic_init2(a, N);
        qadic_init2(b, N);
        qadic_init2(s, N);
        qadic_init2(s1, N);
        qadic_init2(s2, N);
        qadic_init2(lhs, N);
        qadic_init2(rhs, N);

        qadic_randtest_int(a, state, ctx);
        qadic_randtest_int(b, state, ctx);
        e = n_randint(state, 10) % d;

        qadic_mul(s, a, b, ctx);
        qadic_frobenius(lhs, s, e, ctx);
        qadic_frobenius(s1, a, e, ctx);
        qadic_frobenius(s2, b, e, ctx);
        qadic_mul(rhs, s1, s2, ctx);

        result = (qadic_equal(lhs, rhs));
        if (!result)
        {
            flint_printf("FAIL (sigma(a*b) = sigma(a) * sigma(b)):\n\n");
            flint_printf("a = "), qadic_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), qadic_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("s = "), qadic_print_pretty(s, ctx), flint_printf("\n");
            flint_printf("s1 = "), qadic_print_pretty(s1, ctx), flint_printf("\n");
            flint_printf("s2 = "), qadic_print_pretty(s2, ctx), flint_printf("\n");
            flint_printf("lhs = "), qadic_print_pretty(lhs, ctx), flint_printf("\n");
            flint_printf("rhs = "), qadic_print_pretty(rhs, ctx), flint_printf("\n");
            flint_printf("e = %wd\n", e);
            fflush(stdout);
            flint_abort();
        }

        qadic_clear(a);
        qadic_clear(b);
        qadic_clear(s);
        qadic_clear(s1);
        qadic_clear(s2);
        qadic_clear(lhs);
        qadic_clear(rhs);

        fmpz_clear(p);
        qadic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
