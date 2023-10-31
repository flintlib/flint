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

TEST_FUNCTION_START(padic_randtest, state)
{
    int i, result;

    /* Check randtest() */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong lo, hi, N;
        padic_ctx_t ctx;

        padic_t a;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_randtest(a, state, ctx);

        if (N > 0)
        {
            lo = -((N + 9) / 10);
            hi = N;
        }
        else if (N < 0)
        {
            lo = N - ((-N + 9) / 10);
            hi = N;
        }
        else
        {
            lo = -10;
            hi = 0;
        }

        result = padic_is_zero(a) || (lo <= padic_val(a) && padic_val(a) < hi);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("N = %wd\n", N);
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check randtest_not_zero() */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong lo, hi, N;
        padic_ctx_t ctx;

        padic_t a;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_randtest_not_zero(a, state, ctx);

        if (N > 0)
        {
            lo = -((N + 9) / 10);
            hi = N;
        }
        else if (N < 0)
        {
            lo = N - ((-N + 9) / 10);
            hi = N;
        }
        else
        {
            lo = -10;
            hi = 0;
        }

        result = !padic_is_zero(a) && (lo <= padic_val(a) && padic_val(a) < hi);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("N = %wd\n", N);
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
