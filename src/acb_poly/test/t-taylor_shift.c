/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_taylor_shift, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong prec1, prec2;
        acb_poly_t f, g;
        acb_t c, d, e;

        prec1 = 2 + n_randint(state, 500);
        prec2 = 2 + n_randint(state, 500);

        acb_poly_init(f);
        acb_poly_init(g);

        acb_init(c);
        acb_init(d);
        acb_init(e);

        acb_poly_randtest(f, state, 1 + n_randint(state, 40), 1 + n_randint(state, 500), 10);
        acb_poly_randtest(g, state, 1 + n_randint(state, 20), 1 + n_randint(state, 500), 10);

        if (n_randint(state, 2))
            acb_set_si(c, n_randint(state, 5) - 2);
        else
            acb_randtest(c, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));

        if (n_randint(state, 2))
            acb_set_si(d, n_randint(state, 5) - 2);
        else
            acb_randtest(d, state, 1 + n_randint(state, 500), 1 + n_randint(state, 100));

        acb_add(e, c, d, prec1);

        /* check f(x+c)(x+d) = f(x+c+d) */
        acb_poly_taylor_shift(g, f, e, prec2);
        acb_poly_taylor_shift(f, f, c, prec1);
        acb_poly_taylor_shift(f, f, d, prec1);

        if (!acb_poly_overlaps(f, g))
        {
            flint_printf("FAIL\n\n");

            flint_printf("c = "); acb_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_printd(d, 15); flint_printf("\n\n");

            flint_printf("f = "); acb_poly_printd(f, 15); flint_printf("\n\n");
            flint_printf("g = "); acb_poly_printd(g, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_clear(f);
        acb_poly_clear(g);

        acb_clear(c);
        acb_clear(d);
        acb_clear(e);
    }

    TEST_FUNCTION_END(state);
}
