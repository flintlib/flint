/*
    Copyright (C) 2011, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "padic.h"
#include "padic_mat.h"

int
main(void)
{
    int i, result;

    fmpz_t p;
    slong N;
    padic_ctx_t ctx;
    slong m, n;

    FLINT_TEST_INIT(state);

    flint_printf("sub... ");
    fflush(stdout);

    /* Check aliasing: a = a - b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, b, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(b, m, n, N);
        padic_mat_init2(d, m, n, N);

        padic_mat_randtest(a, state, ctx);
        padic_mat_randtest(b, state, ctx);

        padic_mat_sub(d, a, b, ctx);
        padic_mat_sub(a, a, b, ctx);

        result = (padic_mat_equal(a, d) && padic_mat_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_mat_print(b, ctx), flint_printf("\n");
            flint_printf("d = "), padic_mat_print(d, ctx), flint_printf("\n");
            abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(b);
        padic_mat_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: b = a - b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, b, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(b, m, n, N);
        padic_mat_init2(d, m, n, N);

        padic_mat_randtest(a, state, ctx);
        padic_mat_randtest(b, state, ctx);

        padic_mat_sub(d, a, b, ctx);
        padic_mat_sub(b, a, b, ctx);

        result = (padic_mat_equal(b, d) && padic_mat_is_reduced(b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_mat_print(b, ctx), flint_printf("\n");
            flint_printf("d = "), padic_mat_print(d, ctx), flint_printf("\n");
            abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(b);
        padic_mat_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check aliasing: a = a - a == 0 */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(d, m, n, N);

        padic_mat_randtest(a, state, ctx);

        padic_mat_sub(d, a, a, ctx);
        padic_mat_sub(a, a, a, ctx);

        result = (padic_mat_equal(a, d) && padic_mat_is_zero(a));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("d = "), padic_mat_print(d, ctx), flint_printf("\n");
            abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check commutativity: a - b == -(b - a) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, b, c, d;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(b, m, n, N);
        padic_mat_init2(c, m, n, N);
        padic_mat_init2(d, m, n, N);

        padic_mat_randtest(a, state, ctx);
        padic_mat_randtest(b, state, ctx);

        padic_mat_sub(c, a, b, ctx);
        padic_mat_sub(d, b, a, ctx);
        padic_mat_neg(d, d, ctx);

        result = (padic_mat_equal(c, d) && padic_mat_is_reduced(c, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print_pretty(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_mat_print_pretty(b, ctx), flint_printf("\n");
            flint_printf("c = "), padic_mat_print_pretty(c, ctx), flint_printf("\n");
            flint_printf("d = "), padic_mat_print_pretty(d, ctx), flint_printf("\n");
            abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(b);
        padic_mat_clear(c);
        padic_mat_clear(d);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    /* Check a - 0 == a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_mat_t a, b;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN) 
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        padic_mat_init2(a, m, n, N);
        padic_mat_init2(b, m, n, N);

        padic_mat_randtest(a, state, ctx);

        padic_mat_sub(b, a, b, ctx);

        result = (padic_mat_equal(a, b) && padic_mat_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_mat_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_mat_print(b, ctx), flint_printf("\n");
            abort();
        }

        padic_mat_clear(a);
        padic_mat_clear(b);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

