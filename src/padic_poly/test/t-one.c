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
#include "padic_poly.h"

TEST_FUNCTION_START(padic_poly_one, state)
{
    int i, result;

    padic_ctx_t ctx;
    fmpz_t p;
    slong N;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_poly_t a;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_randtest(a, state, n_randint(state, 100), ctx);
        padic_poly_one(a);

        result = (padic_poly_is_one(a) || N <= 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), padic_poly_print(a, ctx), flint_printf("\n\n");
            flint_printf("N = %wd\n\n", N);
            fflush(stdout);
            flint_abort();
        }

        padic_poly_clear(a);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
