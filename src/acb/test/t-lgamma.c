/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_lgamma, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 500);
        prec2 = prec1 + 30;

        acb_init(a);
        acb_init(b);
        acb_init(c);

        arb_randtest_precise(acb_realref(a), state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
        arb_randtest_precise(acb_imagref(a), state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));

        acb_lgamma(b, a, prec1);

        if (n_randint(state, 4) == 0)
        {
            acb_randtest(c, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
            acb_add(a, a, c, prec1);
            acb_sub(a, a, c, prec1);
        }

        acb_lgamma(c, a, prec2);

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        /* check lgamma(z+1) = lgamma(z) + log(z) */
        acb_log(c, a, prec1);
        acb_add(b, b, c, prec1);

        acb_add_ui(c, a, 1, prec1);
        acb_lgamma(c, c, prec1);

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
