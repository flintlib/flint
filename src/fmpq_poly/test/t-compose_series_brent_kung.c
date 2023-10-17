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
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_compose_series_brent_kung, state)
{
    int i, result;

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;
        slong n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpq_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpq_poly_compose_series_brent_kung(f, g, h, n);
        fmpq_poly_compose_series_brent_kung(g, g, h, n);

        result = (fmpq_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL (aliasing 1):\n");
            fmpq_poly_print(f), flint_printf("\n\n");
            fmpq_poly_print(g), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;
        slong n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpq_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpq_poly_compose_series_brent_kung(f, g, h, n);
        fmpq_poly_compose_series_brent_kung(h, g, h, n);

        result = (fmpq_poly_equal(f, h));
        if (!result)
        {
            flint_printf("FAIL (aliasing 2):\n");
            fmpq_poly_print(f), flint_printf("\n\n");
            fmpq_poly_print(h), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Compare with compose */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h, s, t;
        slong n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpq_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpq_poly_compose(s, g, h);
        fmpq_poly_truncate(s, n);
        fmpq_poly_compose_series_brent_kung(f, g, h, n);

        result = (fmpq_poly_equal(f, s));
        if (!result)
        {
            flint_printf("FAIL (comparison):\n");
            flint_printf("n = %wd\n", n);
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("h = "), fmpq_poly_print(h), flint_printf("\n\n");
            flint_printf("f = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("s = "), fmpq_poly_print(s), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}
