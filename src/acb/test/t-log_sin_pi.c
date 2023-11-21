/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_log_sin_pi, state)
{
    slong iter;

    /* check functional equation S(s+1) = S(s) + log(-s) - log(s) */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, s1, r, r1, t;
        slong prec;

        acb_init(s);
        acb_init(s1);
        acb_init(r);
        acb_init(r1);
        acb_init(t);

        prec = 2 + n_randint(state, 500);
        acb_randtest(s, state, 1 + n_randint(state, 500), 1 + n_randint(state, 10));
        acb_randtest(r, state, 1 + n_randint(state, 500), 1 + n_randint(state, 10));
        acb_randtest(r1, state, 1 + n_randint(state, 500), 1 + n_randint(state, 10));

        acb_log_sin_pi(r, s, prec);

        acb_add_ui(s1, s, 1, prec);

        if (n_randint(state, 4) == 0)
        {
            acb_randtest(t, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
            acb_add(s1, s1, t, prec);
            acb_sub(s1, s1, t, prec);
        }

        acb_log_sin_pi(r1, s1, prec);

        acb_log(t, s, prec);
        acb_sub(r, r, t, prec);
        acb_neg(t, s);
        acb_log(t, t, prec);
        acb_add(r, r, t, prec);

        if (!acb_overlaps(r, r1))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("s  = "); acb_printd(s, 30); flint_printf("\n\n");
            flint_printf("s1 = "); acb_printd(s1, 30); flint_printf("\n\n");
            flint_printf("r  = "); acb_printd(r, 30); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(s);
        acb_clear(s1);
        acb_clear(r);
        acb_clear(r1);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
