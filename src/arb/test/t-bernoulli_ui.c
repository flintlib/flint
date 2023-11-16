/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_bernoulli_ui, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t b1, b2;
        ulong n;
        slong prec1, prec2, acc1, acc2;

        n = n_randint(state, 10000);
        prec1 = 2 + n_randint(state, 10000);
        prec2 = prec1 + 100;

        arb_init(b1);
        arb_init(b2);

        arb_bernoulli_ui(b1, n, prec1);
        arb_bernoulli_ui(b2, n, prec2);

        if (!arb_overlaps(b1, b2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("b1 = "); arb_print(b1); flint_printf("\n\n");
            flint_printf("b2 = "); arb_print(b2); flint_printf("\n\n");
            flint_abort();
        }

        acc1 = arb_rel_accuracy_bits(b1);
        acc2 = arb_rel_accuracy_bits(b2);

        if (acc1 < prec1 - 2 || acc2 < prec2 - 2)
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec1 = %wd\n", prec1);
            flint_printf("prec2 = %wd\n", prec2);
            flint_printf("b1 = "); arb_printd(b1, prec1 / 3.33); flint_printf("\n\n");
            flint_printf("b2 = "); arb_printd(b2, prec2 / 3.33); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(b1);
        arb_clear(b2);
    }

    TEST_FUNCTION_END(state);
}
