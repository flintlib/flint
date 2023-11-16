/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_ci, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
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
                arb_hypgeom_ci(s, x, prec1);
                break;
            case 1:
                _arb_hypgeom_ci_2f3(s, x, n_randint(state, prec1), prec1, prec1);
                break;
            default:
                _arb_hypgeom_ci_asymp(s, x, n_randint(state, prec1 / 2), prec1);
                break;
        }

        switch (n_randint(state, 3))
        {
            case 0:
                arb_hypgeom_ci(t, x, prec2);
                break;
            case 1:
                _arb_hypgeom_ci_2f3(t, x, n_randint(state, prec2), prec2, prec2);
                break;
            default:
                _arb_hypgeom_ci_asymp(t, x, n_randint(state, prec2 / 2), prec2);
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

    TEST_FUNCTION_END(state);
}
