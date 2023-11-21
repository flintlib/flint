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

TEST_FUNCTION_START(acb_exp_invexp, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, d, e;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(e);

        acb_randtest(a, state, 1 + n_randint(state, 200), 3);
        acb_randtest(b, state, 1 + n_randint(state, 200), 3);
        acb_randtest(c, state, 1 + n_randint(state, 200), 3);

        prec = 2 + n_randint(state, 200);

        acb_exp_invexp(b, c, a, prec);

        acb_mul(b, b, c, prec);
        acb_one(c);

        if (!acb_contains(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        acb_clear(e);
    }

    TEST_FUNCTION_END(state);
}
