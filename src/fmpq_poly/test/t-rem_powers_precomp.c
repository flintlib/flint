/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_rem_powers_precomp, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of q and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, r;
        fmpq_poly_powers_precomp_t binv;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(r);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);

        fmpq_poly_powers_precompute(binv, b);

        fmpq_poly_rem_powers_precomp(r, a, b, binv);
        fmpq_poly_rem_powers_precomp(a, a, b, binv);

        result = (fmpq_poly_equal(r, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("r = "), fmpq_poly_debug(r), flint_printf("\n\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_powers_clear(binv);

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(r);
    }

    /* Check aliasing of q and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, r;
        fmpq_poly_powers_precomp_t binv;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(r);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);

        fmpq_poly_powers_precompute(binv, b);

        fmpq_poly_rem_powers_precomp(r, a, b, binv);
        fmpq_poly_rem_powers_precomp(b, a, b, binv);

        result = (fmpq_poly_equal(r, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("r = "), fmpq_poly_debug(r), flint_printf("\n\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_powers_clear(binv);

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(r);
    }

    /* Compare with divrem */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, q, r2, r;
        fmpq_poly_powers_precomp_t binv;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);
        fmpq_poly_init(r2);
        fmpq_poly_init(r);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);

        fmpq_poly_powers_precompute(binv, b);

        fmpq_poly_divrem(q, r, a, b);
        fmpq_poly_rem_powers_precomp(r2, a, b, binv);
        fmpq_poly_canonicalise(r2);

        result = (fmpq_poly_equal(r, r2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a  = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b  = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("q  = "), fmpq_poly_debug(q), flint_printf("\n\n");
            flint_printf("r  = "), fmpq_poly_debug(r), flint_printf("\n\n");
            flint_printf("r2 = "), fmpq_poly_debug(r2), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_powers_clear(binv);

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
        fmpq_poly_clear(r2);
        fmpq_poly_clear(r);
    }

    TEST_FUNCTION_END(state);
}
