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
#include "fmpz_poly.h"
#include "padic_poly.h"

TEST_FUNCTION_START(padic_poly_init_realloc_clear, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_poly_t a;

        padic_poly_init2(a, n_randint(state, 100), PADIC_DEFAULT_PREC);
        padic_poly_clear(a);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong N;
        fmpz_t p;
        padic_poly_t a;

        fmpz_init_set_ui(p, 7);

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;

        padic_poly_init2(a, n_randint(state, 100), N);
        padic_poly_realloc(a, n_randint(state, 100), p);
        padic_poly_clear(a);

        fmpz_clear(p);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        padic_ctx_t ctx;
        fmpz_t p;
        slong N;

        padic_poly_t a;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));
        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;
        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_poly_init2(a, 0, N);
        padic_poly_randtest(a, state, n_randint(state, 100), ctx);
        padic_poly_clear(a);

        padic_ctx_clear(ctx);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
