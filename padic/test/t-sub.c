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

    flint_printf("sub... ");
    fflush(stdout);

    

    /* Check aliasing: a = a - b */
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

        padic_sub(d, a, b, ctx);
        padic_sub(a, a, b, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL (alias a = a - b):\n\n");
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

    /* Check aliasing: b = a - b */
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

        padic_sub(d, a, b, ctx);
        padic_sub(b, a, b, ctx);

        result = (padic_equal(b, d));
        if (!result)
        {
            flint_printf("FAIL (alias b = a - b):\n\n");
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

    /* Check aliasing: a = a - a */
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

        padic_sub(d, a, a, ctx);
        padic_sub(a, a, a, ctx);

        result = (padic_equal(a, d));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that a - b == -(b - a) */
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

        padic_sub(c, a, b, ctx);
        padic_sub(d, b, a, ctx);
        padic_neg(d, d, ctx);

        result = (padic_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("c = "), padic_print(c, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that (a - b) - c == a - (b + c) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b, c, d, e;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        padic_init2(c, N);
        padic_init2(d, N);
        padic_init2(e, N);

        padic_randtest(a, state, ctx);
        padic_randtest(b, state, ctx);
        padic_randtest(c, state, ctx);

        padic_sub(d, a, b, ctx);
        padic_sub(d, d, c, ctx);

        padic_add(e, b, c, ctx);
        padic_sub(e, a, e, ctx);

        result = (padic_equal(d, e));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("c = "), padic_print(c, ctx), flint_printf("\n");
            flint_printf("d = "), padic_print(d, ctx), flint_printf("\n");
            flint_printf("e = "), padic_print(e, ctx), flint_printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);
        padic_clear(c);
        padic_clear(d);
        padic_clear(e);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check that a - 0 == a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);

        padic_randtest(a, state, ctx);
        padic_zero(b);

        padic_sub(b, a, b, ctx);

        result = (padic_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            abort();
        }

        padic_clear(a);
        padic_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

