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

TEST_FUNCTION_START(padic_mul, state)
{
    int i, result;

    /* Check aliasing: a = a * b (mod p^N) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);

        padic_mul(d, a, b, ctx);
        padic_mul(a, a, b, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (alias a = a*b):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: b = a * b (mod p^N) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);

        padic_mul(d, a, b, ctx);
        padic_mul(b, a, b, ctx);

        result = (padic_equal(b, d));
        if (!result)
        {
            flint_printf("FAIL (alias b = a*b):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: a = a * a (mod p^N) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);

        padic_mul(d, a, a, ctx);
        padic_mul(a, a, a, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (alias a = a*a):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that a * b == b * a (mod p^N) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, c, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);

        padic_mul(c, a, b, ctx);
        padic_mul(d, b, a, ctx);

        result = (padic_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL (a*b = b*a):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("c = "), padic_print(c, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that (a * b) * c == a * (b * c) (mod p^N) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, c, lhs1, lhs2, rhs1, rhs2;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);
        padic_randtest(c, state, ctx);

        padic_init2(lhs1, N - padic_val(c));
        padic_init2(lhs2, N);
        padic_init2(rhs1, N - padic_val(a));
        padic_init2(rhs2, N);

        padic_mul(lhs1, a, b, ctx);
        padic_mul(lhs2, lhs1, c, ctx);
        padic_mul(rhs1, b, c, ctx);
        padic_mul(rhs2, a, rhs1, ctx);

        result = (padic_equal(lhs2, rhs2));
        if (!result)
        {
            flint_printf("FAIL ((a*b)*c = a*(b*c) mod p^N):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("c = "), padic_print(c, ctx), flint_printf("\n");
            flint_printf("lhs1 = "), padic_print(lhs1, ctx), flint_printf("\n");
            flint_printf("lhs2 = "), padic_print(lhs2, ctx), flint_printf("\n");
            flint_printf("rhs1 = "), padic_print(rhs1, ctx), flint_printf("\n");
            flint_printf("rhs2 = "), padic_print(rhs2, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(lhs1);
        padic_clear(lhs2);
        padic_clear(rhs1);
        padic_clear(rhs2);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that a * 1 == a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, c;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, 1);
        padic_init2(c, N);

        padic_randtest(a, state, ctx);

        padic_one(b);
        padic_mul(c, a, b, ctx);

        result = (padic_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (a*1 = a):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("c = "), padic_print(c, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
