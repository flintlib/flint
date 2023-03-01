/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("add_si....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_poly_t a, b, c, d;
        slong v;

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        arb_poly_randtest(a, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);
        arb_poly_randtest(b, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);
        arb_poly_randtest(c, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);
        arb_poly_randtest(d, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);

        v = n_randtest(state);

        arb_poly_set_si(b, v);

        arb_poly_add(c, a, b, 2 + n_randint(state, 200));
        arb_poly_add_si(d, a, v, 2 + n_randint(state, 200));

        if (!arb_poly_overlaps(c, d))
        {
            flint_printf("FAIL\n\n");
            flint_abort();
        }

        arb_poly_add_si(a, a, v, 2 + n_randint(state, 200));

        if (!arb_poly_overlaps(a, d))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

