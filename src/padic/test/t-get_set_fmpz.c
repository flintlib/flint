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

TEST_FUNCTION_START(padic_get_set_fmpz, state)
{
    int i, result;

    /* Check that Zp(ZZ(x)) == x. */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong N;
        padic_ctx_t ctx;

        padic_t a, b;
        fmpz_t c;

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        N = n_randint(state, PADIC_TEST_PREC_MAX - PADIC_TEST_PREC_MIN)
            + PADIC_TEST_PREC_MIN;

        padic_ctx_init(ctx, p, FLINT_MAX(0, N-10), FLINT_MAX(0, N+10), PADIC_SERIES);

        padic_init2(a, N);
        padic_init2(b, N);
        fmpz_init(c);

        padic_randtest_int(a, state, ctx);

        padic_get_fmpz(c, a, ctx);
        padic_set_fmpz(b, c, ctx);

        result = (padic_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), padic_print(a, ctx), flint_printf("\n");
            flint_printf("b = "), padic_print(b, ctx), flint_printf("\n");
            flint_printf("c = "), fmpz_print(c), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(b);
        fmpz_clear(c);

        fmpz_clear(p);
        padic_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
