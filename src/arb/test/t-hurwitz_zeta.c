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
#include "acb.h"

TEST_FUNCTION_START(arb_hurwitz_zeta, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        acb_t d, e, f;
        slong prec;

        prec = 2 + n_randint(state, 300);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        acb_init(d);
        acb_init(e);
        acb_init(f);

        arb_randtest_precise(a, state, 1 + n_randint(state, 300), 5);
        arb_randtest_precise(b, state, 1 + n_randint(state, 300), 5);
        arb_randtest_precise(c, state, 1 + n_randint(state, 300), 5);

        acb_set_arb(d, a);
        acb_set_arb(e, b);

        arb_hurwitz_zeta(c, a, b, prec);
        acb_hurwitz_zeta(f, d, e, prec);

        if (!arb_overlaps(c, acb_realref(f)) ||
            (arb_is_finite(c) && !arb_contains_zero(acb_imagref(f))))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_printd(d, 15); flint_printf("\n\n");
            flint_printf("e = "); acb_printd(e, 15); flint_printf("\n\n");
            flint_printf("f = "); acb_printd(f, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        acb_clear(d);
        acb_clear(e);
        acb_clear(f);
    }

    TEST_FUNCTION_END(state);
}
