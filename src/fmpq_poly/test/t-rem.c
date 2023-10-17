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

TEST_FUNCTION_START(fmpq_poly_rem, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of r and a */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, r;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(r);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 200);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 200);

        fmpq_poly_rem(r, a, b);
        fmpq_poly_rem(a, a, b);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(r) ? 0 : 2;
        result = (fmpq_poly_equal(r, a) && !cflags);
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

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(r);
    }

    /* Check aliasing of r and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, r;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(r);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 200);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 200);

        fmpq_poly_rem(r, a, b);
        fmpq_poly_rem(b, a, b);

        cflags |= fmpq_poly_is_canonical(b) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(r) ? 0 : 2;
        result = (fmpq_poly_equal(r, b) && !cflags);
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

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(r);
    }

    /* Compare with divrem */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, q, r, r2;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(q);
        fmpq_poly_init(r);
        fmpq_poly_init(r2);
        fmpq_poly_randtest(a, state, n_randint(state, 50), 200);
        fmpq_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 200);

        fmpq_poly_divrem(q, r, a, b);
        fmpq_poly_rem(r2, a, b);

        cflags |= fmpq_poly_is_canonical(q)  ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(r2) ? 0 : 2;
        result = (fmpq_poly_equal(r, r2) && !cflags);
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

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(q);
        fmpq_poly_clear(r);
        fmpq_poly_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
