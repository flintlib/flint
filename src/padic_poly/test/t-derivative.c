/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Sebastian Pancratz

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

TEST_FUNCTION_START(padic_poly_derivative, state)
{
    int i, result;

    padic_ctx_t ctx;
    fmpz_t p;
    slong N;

    /* Aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_poly_t a, b, c;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_init2(b, 0, N);
        padic_poly_init2(c, 0, N);

        padic_poly_randtest(a, state, n_randint(state, 100), ctx);
        padic_poly_set(b, a, ctx);

        padic_poly_derivative(c, b, ctx);
        padic_poly_derivative(b, b, ctx);

        result = (padic_poly_equal(b, c) && padic_poly_is_reduced(b, ctx));
        if (!result)
        {
            flint_printf("FAIL (alias):\n");
            flint_printf("a = "), padic_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), padic_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c = "), padic_poly_print(c, ctx), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    /* Compare with derivative over QQ */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_poly_t a, b, c;
        fmpq_poly_t aQQ, bQQ;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_init2(b, 0, N);
        padic_poly_init2(c, 0, N);

        padic_poly_randtest(a, state, n_randint(state, 100), ctx);

        fmpq_poly_init(aQQ);
        fmpq_poly_init(bQQ);

        padic_poly_derivative(b, a, ctx);

        padic_poly_get_fmpq_poly(aQQ, a, ctx);
        fmpq_poly_derivative(bQQ, aQQ);
        padic_poly_set_fmpq_poly(c, bQQ, ctx);

        result = (padic_poly_equal(b, c) && padic_poly_is_reduced(b, ctx));
        if (!result)
        {
            flint_printf("FAIL (cmp with QQ):\n");
            flint_printf("a = "), padic_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("b = "), padic_poly_print(b, ctx), flint_printf("\n\n");
            flint_printf("c = "), padic_poly_print(c, ctx), flint_printf("\n\n");
            flint_printf("aQQ = "), fmpq_poly_print(aQQ), flint_printf("\n\n");
            flint_printf("bQQ = "), fmpq_poly_print(bQQ), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(a);
        padic_poly_clear(b);
        padic_poly_clear(c);
        fmpq_poly_clear(aQQ);
        fmpq_poly_clear(bQQ);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
