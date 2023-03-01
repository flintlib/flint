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

    flint_printf("sub_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arb_poly_t a, b, c, d;
        slong len, prec;

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        arb_poly_randtest(a, state, 1 + n_randint(state, 10), 100, 10);
        arb_poly_randtest(b, state, 1 + n_randint(state, 10), 100, 10);
        arb_poly_randtest(c, state, 1 + n_randint(state, 10), 100, 10);
        arb_poly_randtest(d, state, 1 + n_randint(state, 10), 100, 10);
        prec = 2 + n_randint(state, 100);
        len = n_randint(state, 10);

        arb_poly_sub_series(c, a, b, len, prec);
        arb_poly_sub(d, a, b, prec);
        arb_poly_truncate(d, len);

        if (!arb_poly_equal(c, d))
        {
            flint_printf("FAIL\n\n");
            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_set(d, a);
        arb_poly_sub_series(d, d, b, len, prec);
        if (!arb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 1)\n\n");
            flint_abort();
        }

        arb_poly_set(d, b);
        arb_poly_sub_series(d, a, d, len, prec);
        if (!arb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 2)\n\n");
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

