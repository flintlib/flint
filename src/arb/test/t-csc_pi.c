/*
    Copyright (C) 2017 Fredrik Johansson
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_csc_pi, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, a, b;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 200);
        prec2 = prec1 + 30;

        arb_init(x);
        arb_init(a);
        arb_init(b);

        arb_randtest_special(x, state, 1 + n_randint(state, 300), 100);
        arb_randtest_special(a, state, 1 + n_randint(state, 300), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 300), 100);

        if (n_randint(state, 2))
        {
            arb_csc_pi(a, x, prec1);
        }
        else
        {
            arb_set(a, x);
            arb_csc_pi(a, a, prec1);
        }

        arb_sin_pi(b, x, prec2);
        arb_inv(b, b, prec2);

        /* check consistency */
        if (!arb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); arb_printn(x, 20, 0); flint_printf("\n\n");
            flint_printf("a = "); arb_printn(a, 20, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 20, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(a);
        arb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
