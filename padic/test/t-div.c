/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "long_extras.h"
#include "padic.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("div... ");
    fflush(stdout);

    

    /* Check aliasing: a = a / b (mod p^N) */
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
        padic_randtest_not_zero(b, state, ctx);

        padic_div(d, a, b, ctx);
        padic_div(a, a, b, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (aliasing a = a/b) :\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: b = a / b (mod p^N) */
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
        padic_randtest_not_zero(b, state, ctx);

        padic_div(d, a, b, ctx);
        padic_div(b, a, b, ctx);

        result = (padic_equal(b, d));
        if (!result)
        {
            flint_printf("FAIL (aliasing b = a/b):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: a = a / a (mod p^N) */
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

        padic_randtest_not_zero(a, state, ctx);

        padic_div(d, a, a, ctx);
        padic_div(a, a, a, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (aliasing a = a/a):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that a / 1 == a */
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
        padic_init2(c, N);
        padic_init2(b, FLINT_MAX(N, 1));

        padic_randtest(a, state, ctx);
        padic_one(b);

        padic_div(c, a, b, ctx);

        result = (padic_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (a/1 == a):\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("c = "), padic_print(b, ctx), flint_printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

