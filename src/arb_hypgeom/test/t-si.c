/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("si....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t x, s, t;
        slong prec1, prec2;

        arb_init(x);
        arb_init(s);
        arb_init(t);

        if (n_randint(state, 10) == 0)
        {
            prec1 = 2 + n_randint(state, 2000);
            prec2 = 2 + n_randint(state, 2000);
        }
        else
        {
            prec1 = 2 + n_randint(state, 200);
            prec2 = 2 + n_randint(state, 200);
        }

        arb_randtest(x, state, 2 + n_randint(state, prec1), 2 + n_randint(state, 100));
        arb_randtest(s, state, 2 + n_randint(state, prec1), 2 + n_randint(state, 100));
        arb_randtest(t, state, 2 + n_randint(state, prec1), 2 + n_randint(state, 100));

        switch (n_randint(state, 3))
        {
            case 0:
                arb_hypgeom_si(s, x, prec1);
                break;
            case 1:
                _arb_hypgeom_si_1f2(s, x, n_randint(state, prec1), prec1, prec1);
                break;
            default:
                _arb_hypgeom_si_asymp(s, x, n_randint(state, prec1 / 2), prec1);
                break;
        }

        switch (n_randint(state, 3))
        {
            case 0:
                arb_hypgeom_si(t, x, prec2);
                break;
            case 1:
                _arb_hypgeom_si_1f2(t, x, n_randint(state, prec2), prec2, prec2);
                break;
            default:
                _arb_hypgeom_si_asymp(t, x, n_randint(state, prec2 / 2), prec2);
                break;
        }

        if (!arb_overlaps(s, t))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("s = "); arb_printn(s, 100, 0); flint_printf("\n\n");
            flint_printf("t = "); arb_printn(t, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(s);
        arb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
