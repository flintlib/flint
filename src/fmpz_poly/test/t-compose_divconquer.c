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
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_compose_divconquer, state)
{
    int i, result;

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpz_poly_compose_divconquer(f, g, h);
        fmpz_poly_compose_divconquer(g, g, h);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(g), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpz_poly_compose_divconquer(f, g, h);
        fmpz_poly_compose_divconquer(h, g, h);

        result = (fmpz_poly_equal(f, h));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(h), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    /* Compare with the default method */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f1, f2, g, h;

        fmpz_poly_init(f1);
        fmpz_poly_init(f2);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);

        fmpz_poly_compose_divconquer(f1, g, h);
        fmpz_poly_compose(f2, g, h);

        result = (fmpz_poly_equal(f1, f2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(f1), flint_printf("\n\n");
            fmpz_poly_print(f2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f1);
        fmpz_poly_clear(f2);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    TEST_FUNCTION_END(state);
}
