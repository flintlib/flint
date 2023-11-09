/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "padic.h"

TEST_FUNCTION_START(padic_inv, state)
{
    int i, result;

/* PRIME p = 2 ***************************************************************/

    /* Check aliasing: a = a^{-1} (mod p^N) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
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

        padic_randtest_not_zero(a, state, ctx);

        padic_inv(d, a, ctx);
        padic_inv(a, a, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
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

    /* Check that correct only mod p^{N} */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, d;
        slong v;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;

        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(d, N);

        padic_randtest_not_zero(a, state, ctx);
        v = padic_val(a);

        if (-v < N) /* Otherwise, no precision left */
        {
            slong N2 = N - FLINT_ABS(v);

            padic_prec(d) = N2;

            padic_inv(b, a, ctx);
            padic_mul(d, a, b, ctx);

            result = (padic_is_one(d));
            if (!result)
            {
                flint_printf("FAIL (a * a^{-1} == 1):\n\n");
                flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
                flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
                flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

/* PRIME p > 2 ***************************************************************/

    /* Check aliasing: a = a^{-1} (mod p^N) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
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

        padic_randtest_not_zero(a, state, ctx);

        padic_inv(d, a, ctx);
        padic_inv(a, a, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
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

    /* Check that correct only mod p^{N} */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, d;
        slong v;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;

        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(d, N);

        padic_randtest_not_zero(a, state, ctx);
        v = padic_val(a);

        if (-v < N) /* Otherwise, no precision left */
        {
            slong N2 = N - FLINT_ABS(v);

            padic_prec(d) = N2;

            padic_inv(b, a, ctx);
            padic_mul(d, a, b, ctx);

            result = (padic_is_one(d));
            if (!result)
            {
                flint_printf("FAIL (a * a^{-1} == 1):\n\n");
                flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
                flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
                flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
