/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_sinc, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, d;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);

        acb_randtest(a, state, 1 + n_randint(state, 400), 2 + n_randint(state, 100));
        acb_randtest(b, state, 1 + n_randint(state, 400), 2 + n_randint(state, 100));
        acb_randtest(c, state, 1 + n_randint(state, 400), 2 + n_randint(state, 100));

        acb_sinc(b, a, 2 + n_randint(state, 400));
        acb_sin(c, a, 2 + n_randint(state, 400));
        acb_mul(d, b, a, 2 + n_randint(state, 400));

        /* check consistency */
        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
            flint_printf("d = "); acb_printd(d, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_set(c, a);
        acb_sinc(c, c, 2 + n_randint(state, 400));

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(c, 30); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
