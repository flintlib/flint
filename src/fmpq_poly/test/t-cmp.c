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

TEST_FUNCTION_START(fmpq_poly_cmp, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f;

        fmpq_poly_init(f);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);

        result = (fmpq_poly_cmp(f, f) == 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
    }

    /*
       Check transitivity, i.e. f <= g <= h implies f <= h, that is
       NOT (f <= g <= h) OR f <= h
     */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(g, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(h, state, n_randint(state, 100), 200);

        result = !(fmpq_poly_cmp(f, g) <= 0) || !(fmpq_poly_cmp(g, h) <= 0)
            || (fmpq_poly_cmp(f, h) <= 0);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n");
            fmpq_poly_debug(g), flint_printf("\n");
            fmpq_poly_debug(h), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Check that <, ==, or > */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(g, state, n_randint(state, 100), 200);

        result = (fmpq_poly_cmp(f, g) < 0) || (fmpq_poly_equal(f, g))
            || (fmpq_poly_cmp(f, g) > 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(f), flint_printf("\n");
            fmpq_poly_debug(g), flint_printf("\n");
            flint_printf("cmp(f,g) = %d\n", fmpq_poly_cmp(f, g));
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
