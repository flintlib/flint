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

TEST_FUNCTION_START(padic_teichmuller, state)
{
    int i, result;

    /* Check aliasing (x 1,000) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, c;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);

        padic_randtest_int(a, state, ctx);
        padic_set(b, a, ctx);

        padic_teichmuller(c, b, ctx);
        padic_teichmuller(b, b, ctx);

        result = (padic_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
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

    /* Check x^p == x for word-sized p (x 10,000)*/
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong prime, N;
        padic_ctx_t ctx;

        padic_t a, b, c;

        prime = n_randprime(state, 2 + n_randint(state, SMALL_FMPZ_BITCOUNT_MAX), 0);
        fmpz_init_set_ui(p, prime);
        N = n_randint(state, PADIC_TEST_PREC_MAX);
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);

        padic_randtest_int(a, state, ctx);

        padic_teichmuller(b, a, ctx);
        padic_pow_si(c, b, fmpz_get_si(p), ctx);

        result = (padic_equal(b, c));
        if (!result)
        {
            flint_printf("FAIL (x^p == x):\n\n");
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
