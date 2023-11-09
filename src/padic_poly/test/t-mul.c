/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq_poly.h"
#include "padic_poly.h"

TEST_FUNCTION_START(padic_poly_mul, state)
{
    int i, result;

    padic_ctx_t ctx;
    fmpz_t p;
    slong N;

    /* Check aliasing of a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_poly_t a, b, c;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_init2(b, 0, N);
        padic_poly_init2(c, 0, N);

        padic_poly_randtest(b, state, n_randint(state, 50), ctx);
        padic_poly_randtest(c, state, n_randint(state, 50), ctx);

        padic_poly_mul(a, b, c, ctx);
        padic_poly_mul(b, b, c, ctx);

        result = (padic_poly_equal(a, b) && padic_poly_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and b):\n");
            padic_poly_print(a, ctx), flint_printf("\n\n");
            padic_poly_print(b, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_poly_t a, b, c;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_init2(b, 0, N);
        padic_poly_init2(c, 0, N);

        padic_poly_randtest(b, state, n_randint(state, 50), ctx);
        padic_poly_randtest(c, state, n_randint(state, 50), ctx);

        padic_poly_mul(a, b, c, ctx);
        padic_poly_mul(c, b, c, ctx);

        result = (padic_poly_equal(a, c) && padic_poly_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and c):\n");
            padic_poly_print(a, ctx), flint_printf("\n\n");
            padic_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* Check (b * c) + (b * d) = b * (c + d) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        padic_poly_t a1, a2, b, c, d, t;
        slong v;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(b, 0, N);
        padic_poly_init2(c, 0, N);
        padic_poly_init2(d, 0, N);
        padic_poly_init2(t, 0, N);

        padic_poly_randtest(b, state, n_randint(state, 100), ctx);
        padic_poly_randtest(c, state, n_randint(state, 100), ctx);
        padic_poly_randtest(d, state, n_randint(state, 100), ctx);

        v = FLINT_MIN(b->val, c->val);
        v = FLINT_MIN(v, d->val);
        v = FLINT_MIN(v, 0);

        if (v >= 0 || -v < N)  /* Otherwise, no precision left */
        {
            slong N2 = (v >= 0) ? N : N + v;

            padic_poly_init2(a1, 0, N2);
            padic_poly_init2(a2, 0, N2);

            padic_poly_mul(a1, b, c, ctx);
            padic_poly_mul(t, b, d, ctx);
            padic_poly_add(a1, a1, t, ctx);     /* Lower precision */

            padic_poly_add(t, c, d, ctx);
            padic_poly_mul(a2, b, t, ctx);      /* Lower precision */

            result = (padic_poly_equal(a1, a2) && padic_poly_is_reduced(a1, ctx));
            if (!result)
            {
                flint_printf("FAIL (distributivity):\n");
                flint_printf("p = "), fmpz_print(ctx->p), flint_printf("\n\n");
                flint_printf("N = %wd\n\n", N);
                flint_printf("b = "), padic_poly_print(b, ctx), flint_printf("\n\n");
                flint_printf("c = "), padic_poly_print(c, ctx), flint_printf("\n\n");
                flint_printf("d = "), padic_poly_print(d, ctx), flint_printf("\n\n");
                flint_printf("a1 = "), padic_poly_print(a1, ctx), flint_printf("\n\n");
                flint_printf("a2 = "), padic_poly_print(a2, ctx), flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }

            padic_poly_clear(a1);
            padic_poly_clear(a2);
        }

        padic_poly_clear(b);
        padic_poly_clear(c);
        padic_poly_clear(d);
        padic_poly_clear(t);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* Compare with Q */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_poly_t a, b, c, d;
        fmpq_poly_t x, y, z;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_init2(b, 0, N);
        padic_poly_init2(c, 0, N);
        padic_poly_init2(d, 0, N);

        fmpq_poly_init(x);
        fmpq_poly_init(y);
        fmpq_poly_init(z);

        padic_poly_randtest(b, state, n_randint(state, 50), ctx);
        padic_poly_randtest(c, state, n_randint(state, 50), ctx);

        padic_poly_mul(a, b, c, ctx);

        padic_poly_get_fmpq_poly(y, b, ctx);
        padic_poly_get_fmpq_poly(z, c, ctx);

        fmpq_poly_mul(x, y, z);
        padic_poly_set_fmpq_poly(d, x, ctx);

        result = (padic_poly_equal(a, d) && padic_poly_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL (cmp with Q):\n");
            flint_printf("N = %wd, val(b) = %wd, val(c) = %wd\n", N, b->val, c->val);
            padic_poly_print(c, ctx), flint_printf("\n\n");
            padic_poly_print(d, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);
        padic_poly_clear(d);

        fmpq_poly_clear(x);
        fmpq_poly_clear(y);
        fmpq_poly_clear(z);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
