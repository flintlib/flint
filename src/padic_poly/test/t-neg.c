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
#include "padic_poly.h"

TEST_FUNCTION_START(padic_poly_neg, state)
{
    int i, result;

    fmpz_t p;
    slong N;
    padic_ctx_t ctx;

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

        padic_poly_neg(c, b, ctx);
        padic_poly_neg(b, b, ctx);

        result = (padic_poly_equal(b, c) && padic_poly_is_reduced(b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            padic_poly_print(a, ctx), flint_printf("\n\n");
            padic_poly_print(b, ctx), flint_printf("\n\n");
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

    /* Check that --a = a */
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

        padic_poly_neg(b, a, ctx);
        padic_poly_neg(c, b, ctx);

        result = (padic_poly_equal(a, c) && padic_poly_is_reduced(a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            padic_poly_print(a, ctx), flint_printf("\n\n");
            padic_poly_print(b, ctx), flint_printf("\n\n");
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

    TEST_FUNCTION_END(state);
}
