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

TEST_FUNCTION_START(acb_barnes_g, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, d;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 200);
        prec2 = prec1 + 30;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);

        arb_randtest_precise(acb_realref(a), state, 1 + n_randint(state, 300), 1 + n_randint(state, 8));
        arb_randtest_precise(acb_imagref(a), state, 1 + n_randint(state, 300), 1 + n_randint(state, 8));

        acb_log_barnes_g(b, a, prec1);

        if (n_randint(state, 4) == 0)
        {
            acb_randtest(c, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
            acb_add(a, a, c, prec1);
            acb_sub(a, a, c, prec1);
        }

        acb_log_barnes_g(c, a, prec2);

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        acb_barnes_g(d, a, prec1);
        acb_exp(c, b, prec1);

        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: overlap 2\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        /* check lG(z+1) = lG(z) + lgamma(z) */
        acb_lgamma(c, a, prec1);
        acb_add(b, b, c, prec1);

        acb_add_ui(c, a, 1, prec1);
        acb_log_barnes_g(c, c, prec1);

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
        acb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
