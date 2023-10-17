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

TEST_FUNCTION_START(fmpq_poly_make_monic, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);

        fmpq_poly_make_monic(g, f);
        fmpq_poly_make_monic(f, f);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n\n");
            fmpq_poly_debug(g), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    /* Check that the result of "monic" has "is_monic" return 1 */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f;

        fmpq_poly_init(f);
        fmpq_poly_randtest_not_zero(f, state, n_randint(state, 100) + 1, 200);

        fmpq_poly_make_monic(f, f);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        result = (fmpq_poly_is_monic(f) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
