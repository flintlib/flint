/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_log_hypot, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, r, s;
        slong prec1, prec2, acc1, acc2;

        prec1 = 2 + n_randint(state, 400);
        prec2 = 2 + n_randint(state, 400);

        arb_init(a);
        arb_init(b);
        arb_init(r);
        arb_init(s);

        arb_randtest_special(a, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        arb_randtest_special(r, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));
        arb_randtest_special(s, state, 1 + n_randint(state, 1000), 2 + n_randint(state, 100));

        if (n_randint(state, 2))
        {
            arb_log_hypot(r, a, b, prec1);
        }
        else if (n_randint(state, 2))
        {
            arb_set(r, a);
            arb_log_hypot(r, r, b, prec1);
        }
        else
        {
            arb_set(r, b);
            arb_log_hypot(r, a, r, prec1);
        }

        arb_hypot(s, a, b, prec2);
        arb_log(s, s, prec2);

        /* check consistency */
        if (!arb_overlaps(r, s))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("r = "); arb_printn(r, 50, 0); flint_printf("\n\n");
            flint_printf("s = "); arb_printn(s, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (!arb_is_finite(r) && arb_is_finite(a) && arb_is_finite(b) &&
            (!arb_contains_zero(a) || !arb_contains_zero(b)))
        {
            flint_printf("FAIL: not finite\n\n");
            flint_printf("prec1 = %wd\n\n", prec1);
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("r = "); arb_printn(r, 50, 0); flint_printf("\n\n");
            flint_printf("s = "); arb_printn(s, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        acc1 = arb_rel_accuracy_bits(r);
        acc2 = arb_rel_accuracy_bits(s);

        if (prec2 <= prec1 && acc1 > 0 && acc2 > 0 && acc1 < acc2 - 2)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("prec1 = %wd, acc1 = %wd, acc2 = %wd\n\n", prec1, acc1, acc2);
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("r = "); arb_printn(r, 50, 0); flint_printf("\n\n");
            flint_printf("s = "); arb_printn(s, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (!arf_is_nan(arb_midref(a)) && !arf_is_nan(arb_midref(b)))
        {
            if ((!arb_is_finite(a) && !arb_contains_zero(a)) ||
                (!arb_is_finite(b) && !arb_contains_zero(b)))
            {
                if (arf_is_nan(arb_midref(r)) || arb_is_finite(r) ||
                    !arb_is_positive(r))
                {
                    flint_printf("FAIL: infinite\n\n");
                    flint_printf("prec1 = %wd, acc1 = %wd, acc2 = %wd\n\n", prec1, acc1, acc2);
                    flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n\n");
                    flint_printf("b = "); arb_printn(b, 50, 0); flint_printf("\n\n");
                    flint_printf("r = "); arb_printn(r, 50, 0); flint_printf("\n\n");
                    flint_printf("s = "); arb_printn(s, 50, 0); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(r);
        arb_clear(s);
    }

    TEST_FUNCTION_END(state);
}
