/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_erfc, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        slong prec1, prec2, prec3, prec4;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);
        prec3 = 2 + n_randint(state, 1000);
        prec4 = 2 + n_randint(state, 1000);

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest_special(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest_special(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest_special(c, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        switch (n_randint(state, 4))
        {
            case 0:
                acb_hypgeom_erf_asymp(b, a, 1, prec1, prec3);
                break;
            case 1:
                acb_hypgeom_erf(b, a, prec1);
                acb_sub_ui(b, b, 1, prec1);
                acb_neg(b, b);
                break;
            default:
                acb_hypgeom_erfc(b, a, prec1);
        }

        switch (n_randint(state, 4))
        {
            case 0:
                acb_hypgeom_erf_asymp(c, a, 1, prec2, prec4);
                break;
            case 1:
                acb_hypgeom_erf(c, a, prec2);
                acb_sub_ui(c, c, 1, prec2);
                acb_neg(c, c);
                break;
            default:
                acb_hypgeom_erfc(c, a, prec2);
        }

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
