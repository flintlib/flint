/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_rising_ui_series, state)
{
    slong iter;

    /* check rf(f, a) * rf(f + a, b) = rf(f, a + b) */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong bits, trunc;
        ulong a, b;
        arb_poly_t f, g, h1, h2, h1h2, h3;

        bits = 2 + n_randint(state, 200);
        trunc = 1 + n_randint(state, 20);
        a = n_randint(state, 10);
        b = n_randint(state, 10);

        arb_poly_init(f);
        arb_poly_init(g);
        arb_poly_init(h1);
        arb_poly_init(h2);
        arb_poly_init(h1h2);
        arb_poly_init(h3);

        arb_poly_randtest(f, state, 1 + n_randint(state, 20), bits, 4);
        arb_poly_set(g, f);

        /* g = f + 1 */
        if (g->length == 0)
        {
            arb_poly_fit_length(g, 1);
            arb_set_ui(g->coeffs, a);
            _arb_poly_set_length(g, 1);
            _arb_poly_normalise(g);
        }
        else
        {
            arb_add_ui(g->coeffs, g->coeffs, a, bits);
            _arb_poly_normalise(g);
        }

        arb_poly_rising_ui_series(h1, f, a, trunc, bits);
        arb_poly_rising_ui_series(h2, g, b, trunc, bits);
        arb_poly_rising_ui_series(h3, f, a + b, trunc, bits);

        arb_poly_mullow(h1h2, h1, h2, trunc, bits);

        if (!arb_poly_overlaps(h1h2, h3))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits = %wd\n", bits);
            flint_printf("trunc = %wd\n", trunc);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", a);

            flint_printf("f = "); arb_poly_printd(f, 15); flint_printf("\n\n");
            flint_printf("g = "); arb_poly_printd(g, 15); flint_printf("\n\n");
            flint_printf("h1 = "); arb_poly_printd(h1, 15); flint_printf("\n\n");
            flint_printf("h2 = "); arb_poly_printd(h2, 15); flint_printf("\n\n");
            flint_printf("h1h2 = "); arb_poly_printd(h1h2, 15); flint_printf("\n\n");
            flint_printf("h3 = "); arb_poly_printd(h3, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_rising_ui_series(f, f, a, trunc, bits);

        if (!arb_poly_equal(f, h1))
        {
            flint_printf("FAIL (aliasing)\n\n");

            flint_printf("bits = %wd\n", bits);
            flint_printf("trunc = %wd\n", trunc);
            flint_printf("a = %wu\n", a);
            flint_printf("b = %wu\n", a);

            flint_printf("f = "); arb_poly_printd(f, 15); flint_printf("\n\n");
            flint_printf("h1 = "); arb_poly_printd(h1, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_clear(f);
        arb_poly_clear(g);
        arb_poly_clear(h1);
        arb_poly_clear(h2);
        arb_poly_clear(h1h2);
        arb_poly_clear(h3);
    }

    TEST_FUNCTION_END(state);
}
