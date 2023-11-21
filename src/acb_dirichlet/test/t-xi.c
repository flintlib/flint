/*
    Copyright (C) 2018 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_xi, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest_precise(a, state, 1 + n_randint(state, 1000), 3);

        acb_dirichlet_xi(b, a, prec1);

        if (n_randint(state, 2))
        {
            acb_dirichlet_xi(c, a, prec2);
        }
        else  /* test aliasing */
        {
            acb_set(c, a);
            acb_dirichlet_xi(c, c, prec2);
        }

        if (!acb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        /* check xi(s) = xi(1-s) */
        acb_sub_ui(c, a, 1, prec1);
        acb_neg(c, c);
        acb_dirichlet_xi(c, c, prec1);

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
