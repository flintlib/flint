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
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_rescale, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;
        fmpq_t a;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 100);

        fmpq_init(a);
        fmpq_randtest(a, state, 100);

        fmpq_poly_rescale(g, f, a);
        fmpq_poly_rescale(f, f, a);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n");
            fmpq_poly_debug(f), flint_printf("\n\n");
            fmpq_poly_debug(g), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_clear(a);
    }

    /* Check that rescaling by a and then by 1/a is the identity */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;
        fmpq_t a;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 100);

        fmpq_init(a);
        fmpq_randtest_not_zero(a, state, 100);

        fmpq_poly_rescale(g, f, a);
        fmpq_inv(a, a);
        fmpq_poly_rescale(g, g, a);

        cflags |= fmpq_poly_is_canonical(f) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(g) ? 0 : 2;
        result = (fmpq_poly_equal(f, g) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (composition of a and 1/a):\n");
            fmpq_poly_debug(f), flint_printf("\n\n");
            fmpq_poly_debug(g), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_clear(a);
    }

    TEST_FUNCTION_END(state);
}
