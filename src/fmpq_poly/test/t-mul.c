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

TEST_FUNCTION_START(fmpq_poly_mul, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(b, state, n_randint(state, 50), 500);
        fmpq_poly_randtest(c, state, n_randint(state, 50), 500);

        fmpq_poly_mul(a, b, c);
        fmpq_poly_mul(b, b, c);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(b, state, n_randint(state, 50), 500);
        fmpq_poly_randtest(c, state, n_randint(state, 50), 500);

        fmpq_poly_mul(a, b, c);
        fmpq_poly_mul(c, b, c);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        result = (fmpq_poly_equal(a, c) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    /* Check (b*c)+(b*d) = b*(c+d) */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a1, a2, b, c, d;

        fmpq_poly_init(a1);
        fmpq_poly_init(a2);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_init(d);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 500);
        fmpq_poly_randtest(c, state, n_randint(state, 100), 500);
        fmpq_poly_randtest(d, state, n_randint(state, 100), 500);

        fmpq_poly_mul(a1, b, c);
        fmpq_poly_mul(a2, b, d);
        fmpq_poly_add(a1, a1, a2);

        fmpq_poly_add(c, c, d);
        fmpq_poly_mul(a2, b, c);

        cflags |= fmpq_poly_is_canonical(a1) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(a2) ? 0 : 2;
        result = (fmpq_poly_equal(a1, a2) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a1), flint_printf("\n\n");
            fmpq_poly_debug(a2), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a1);
        fmpq_poly_clear(a2);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpq_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
