/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_rgamma, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(a);
        acb_init(b);
        acb_init(c);

        arb_randtest_precise(acb_realref(a), state, 1 + n_randint(state, 1000), 3);
        arb_randtest_precise(acb_imagref(a), state, 1 + n_randint(state, 1000), 3);

        acb_rgamma(b, a, prec1);
        acb_rgamma(c, a, prec2);

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        /* check 1/gamma(z+1) = 1/gamma(z)/z */
        acb_div(b, b, a, prec1);
        acb_add_ui(c, a, 1, prec1);
        acb_rgamma(c, c, prec1);

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
