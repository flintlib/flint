/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_get_slice, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;
        slong j1, j2;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_set(b, a);

        j1 = n_randint(state, 100);
        j2 = n_randint(state, 100);

        fmpq_poly_get_slice(c, b, j1, j2);
        fmpq_poly_get_slice(b, b, j1, j2);

        cflags |= fmpq_poly_is_canonical(b) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        result = (fmpq_poly_equal(b, c) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    /* Check slice with i >= j is zero */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong j1, j2;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        j2 = n_randint(state, 60);
        j1 = j2 + n_randint(state, 60);

        fmpq_poly_get_slice(b, a, j1, j2);

        cflags |= fmpq_poly_is_canonical(b) ? 0 : 1;
        result = (fmpq_poly_is_zero(b) && !cflags);
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
    }

    /* Check transitivity when j1 <= k1 <= k2 <= j2 */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c, d;
        slong j1, j2, k1, k2;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_init(d);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        j1 = n_randint(state, 20);
        k1 = j1 + n_randint(state, 20);
        k2 = k1 + n_randint(state, 20);
        j2 = k2 + n_randint(state, 20);

        fmpq_poly_get_slice(b, a, j1, j2);
        fmpq_poly_get_slice(c, b, k1, k2);
        fmpq_poly_get_slice(d, a, k1, k2);

        cflags |= fmpq_poly_is_canonical(b) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(d) ? 0 : 4;
        result = (fmpq_poly_equal(c, d) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            fmpq_poly_debug(d), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
        fmpq_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
