/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_central_bin_ui, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t b1, b2, t;
        ulong n;
        slong prec1, prec2, acc1;

        n = n_randtest(state);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(b1);
        arb_init(b2);
        arb_init(t);

        arb_hypgeom_central_bin_ui(b1, n, prec1);

        arb_set_ui(t, n);
        arb_add_ui(t, t, n, prec2);
        arb_add_ui(t, t, 1, prec2);
        arb_gamma(t, t, prec2);
        arb_set_ui(b2, n);
        arb_add_ui(b2, b2, 1, prec2);
        arb_rgamma(b2, b2, prec2);
        arb_mul(b2, b2, b2, prec2);
        arb_mul(b2, b2, t, prec2);

        if (!arb_overlaps(b1, b2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("b1 = "); arb_printn(b1, 50, 0); flint_printf("\n\n");
            flint_printf("b2 = "); arb_printn(b2, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        acc1 = arb_rel_accuracy_bits(b1);

        if (acc1 < prec1 - 2)
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec1 = %wd, acc1 = %wd\n", prec1, acc1);
            flint_printf("b1 = "); arb_printn(b1, prec1 / 3.33, 0); flint_printf("\n\n");
            flint_printf("b2 = "); arb_printn(b2, prec2 / 3.33, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(b1);
        arb_clear(b2);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
