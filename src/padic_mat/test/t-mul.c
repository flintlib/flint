/*
    Copyright (C) 2011, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "padic_mat.h"

TEST_FUNCTION_START(padic_mat_mul, state)
{
    int i, result;

    fmpz_t p;
    slong N;
    padic_ctx_t ctx;
    slong m, n;

    /* Check aliasing: a = a * b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, b, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_VAL_UNIT);

        m = n_randint(state, 10);
        n = m;

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(b, m, n, N);
        padic_mat_init2(d, m, n, N);

        padic_mat_randtest(a, state, ctx);
        padic_mat_randtest(b, state, ctx);

        padic_mat_mul(d, a, b, ctx);
        padic_mat_mul(a, a, b, ctx);

        result = (padic_mat_equal(a, d) && padic_mat_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_mat_print(b, ctx), flint_printf("\n");
            flint_printf("d = "), padic_mat_print(d, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(b);
        padic_mat_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: b = a * b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, b, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_VAL_UNIT);

        m = n_randint(state, 10);
        n = m;

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(b, m, n, N);
        padic_mat_init2(d, m, n, N);

        padic_mat_randtest(a, state, ctx);
        padic_mat_randtest(b, state, ctx);

        padic_mat_mul(d, a, b, ctx);
        padic_mat_mul(b, a, b, ctx);

        result = (padic_mat_equal(b, d) && padic_mat_is_reduced(b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_mat_print(b, ctx), flint_printf("\n");
            flint_printf("d = "), padic_mat_print(d, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(b);
        padic_mat_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: a = a * a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_VAL_UNIT);

        m = n_randint(state, 10);
        n = m;

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(d, m, n, N);

        padic_mat_randtest(a, state, ctx);

        padic_mat_mul(d, a, a, ctx);
        padic_mat_mul(a, a, a, ctx);

        result = (padic_mat_equal(a, d) && padic_mat_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("d = "), padic_mat_print(d, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check identity: a * Id == a, for N > 0 */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, b;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        N = FLINT_MAX(1, N);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_VAL_UNIT);

        m = n_randint(state, 10);
        n = m;

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(b, m, n, N);

        padic_mat_randtest(a, state, ctx);

        padic_mat_one(b);
        padic_mat_mul(b, a, b, ctx);

        result = (padic_mat_equal(a, b) && padic_mat_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL (A * Id == A):\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_mat_print(b, ctx), flint_printf("\n");
            flint_printf("p = "), fmpz_print(p), flint_printf("\n");
            flint_printf("N = %wd\n", N);
            fflush(stdout);
            flint_abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check associativity: (a*b)*c == a*(b*c) mod p^{N-v} */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, b, c, d, e, t1, t2;
        slong k, l, v;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_VAL_UNIT);

        k = n_randint(state, 10);
        l = n_randint(state, 10);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        padic_mat_init2(a, k, l, N);
        padic_mat_init2(b, l, m, N);
        padic_mat_init2(c, m, n, N);

        padic_mat_randtest(a, state, ctx);
        padic_mat_randtest(b, state, ctx);
        padic_mat_randtest(c, state, ctx);

     /* v = min(val(a), val(b), val(c), 0) */
        v = FLINT_MIN(padic_mat_val(a), padic_mat_val(b));
        v = FLINT_MIN(v, padic_mat_val(c));
        v = FLINT_MIN(v, 0);

        if ((v >= 0) || (-v < N))  /* Otherwise, no precision left */
        {
            slong N2 = (v >= 0) ? N : N + v;

            padic_mat_init2(d, k, n, N2);
            padic_mat_init2(e, k, n, N2);
            padic_mat_init2(t1, k, m, N);
            padic_mat_init2(t2, l, n, N);

            padic_mat_mul(t1, a, b, ctx);
            padic_mat_mul(d, t1, c, ctx);
            padic_mat_mul(t2, b, c, ctx);
            padic_mat_mul(e, a, t2, ctx);

            result = (padic_mat_equal(d, e) && padic_mat_is_reduced(d, ctx));
            if (!result)
            {
                flint_printf("FAIL:\n\n");
                flint_printf("a = "), padic_mat_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), padic_mat_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), padic_mat_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), padic_mat_print_pretty(d, ctx), flint_printf("\n");
                flint_printf("e = "), padic_mat_print_pretty(e, ctx), flint_printf("\n");
                flint_printf("t1 = "), padic_mat_print_pretty(t1, ctx), flint_printf("\n");
                flint_printf("t2 = "), padic_mat_print_pretty(t2, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            padic_mat_clear(d);
            padic_mat_clear(e);
            padic_mat_clear(t1);
            padic_mat_clear(t2);
        }

        padic_mat_clear(a);
        padic_mat_clear(b);
        padic_mat_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check distributivity: a(b + c) == ab + ac, precision loss */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, b, c;
        slong l;
        slong v;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_VAL_UNIT);

        l = n_randint(state, 10);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        padic_mat_init2(a, l, m, N);
        padic_mat_init2(b, m, n, N);
        padic_mat_init2(c, m, n, N);

        padic_mat_randtest(a, state, ctx);
        padic_mat_randtest(b, state, ctx);
        padic_mat_randtest(c, state, ctx);

        v = FLINT_MIN(a->val, b->val);
        v = FLINT_MIN(v, c->val);
        v = FLINT_MIN(v, 0);

        if (v >= 0 || -v < N)  /* Otherwise, no precision left */
        {
            slong N2 = (v >= 0) ? N : N + v;

            padic_mat_t lhs, rhs, s, t;

            padic_mat_init2(lhs, l, n, N2);
            padic_mat_init2(rhs, l, n, N2);
            padic_mat_init2(s, m, n, N);
            padic_mat_init2(t, l, n, N2);

            padic_mat_add(s, b, c, ctx);
            padic_mat_mul(lhs, a, s, ctx);

            padic_mat_mul(rhs, a, b, ctx);
            padic_mat_mul(t, a, c, ctx);
            padic_mat_add(rhs, rhs, t, ctx);

            result = (padic_mat_equal(lhs, rhs) && padic_mat_is_reduced(lhs, ctx));
            if (!result)
            {
                flint_printf("FAIL:\n\n");
                flint_printf("Hier...\n");
                flint_printf("l m n = %wd %wd %wd\n", l, m, n);
                flint_printf("N     = %wd\n", N);
                flint_printf("N2    = %wd\n", N2);
                flint_printf("a = "), padic_mat_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), padic_mat_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), padic_mat_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("lhs = "), padic_mat_print_pretty(lhs, ctx), flint_printf("\n");
                flint_printf("rhs = "), padic_mat_print_pretty(rhs, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            padic_mat_clear(lhs);
            padic_mat_clear(rhs);
            padic_mat_clear(s);
            padic_mat_clear(t);
        }

        padic_mat_clear(a);
        padic_mat_clear(b);
        padic_mat_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
