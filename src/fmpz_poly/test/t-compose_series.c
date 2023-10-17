/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

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

TEST_FUNCTION_START(fmpz_poly_compose_series, state)
{
    int i, result;

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h;
        slong n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpz_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpz_poly_compose_series(f, g, h, n);
        fmpz_poly_compose_series(g, g, h, n);

        result = (fmpz_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL (aliasing 1):\n");
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
        slong n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpz_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpz_poly_compose_series(f, g, h, n);
        fmpz_poly_compose_series(h, g, h, n);

        result = (fmpz_poly_equal(f, h));
        if (!result)
        {
            flint_printf("FAIL (aliasing 2):\n");
            fmpz_poly_print(f), flint_printf("\n\n");
            fmpz_poly_print(h), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
    }

    /* Compare with compose */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h, s, t;
        slong n;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(s);
        fmpz_poly_init(t);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpz_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 10);

        fmpz_poly_compose(s, g, h);
        fmpz_poly_truncate(s, n);
        fmpz_poly_compose_series(f, g, h, n);

        result = (fmpz_poly_equal(f, s));
        if (!result)
        {
            flint_printf("FAIL (comparison):\n");
            flint_printf("n = %wd\n", n);
            flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n\n");
            flint_printf("h = "), fmpz_poly_print(h), flint_printf("\n\n");
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("s = "), fmpz_poly_print(s), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(s);
        fmpz_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
