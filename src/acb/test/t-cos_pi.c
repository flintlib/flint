/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_cos_pi, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, y, a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(x);
        acb_init(y);
        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        acb_cos_pi(a, x, prec1);
        acb_cos_pi(b, x, prec2);

        /* check consistency */
        if (!acb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        /* compare with cos */
        arb_const_pi(acb_realref(c), prec1);
        acb_mul_arb(y, x, acb_realref(c), prec1);
        acb_cos(c, y, prec1);

        if (!acb_overlaps(a, c))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_printf("y = "); acb_print(y); flint_printf("\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        acb_cos_pi(x, x, prec1);

        if (!acb_overlaps(a, x))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("x = "); acb_print(x); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
