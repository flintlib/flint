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

TEST_FUNCTION_START(acb_atan, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, a, b;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        acb_init(x);
        acb_init(a);
        acb_init(b);

        acb_randtest_special(x, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        acb_randtest_special(a, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        acb_randtest_special(b, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        acb_atan(a, x, prec1);
        acb_atan(b, x, prec2);

        /* check consistency */
        if (!acb_overlaps(a, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_printd(x, 15); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 15); flint_printf("\n\n");
            flint_abort();
        }

        /* check tan(atan(x)) = x */
        acb_tan(b, b, prec1);

        if (!acb_contains(b, x))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("x = "); acb_printd(x, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_atan(x, x, prec1);

        if (!acb_overlaps(a, x))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n\n");
            flint_printf("x = "); acb_printd(x, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(a);
        acb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
