/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_eta, state)
{
    slong iter;

    for (iter = 0; iter < 250 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z1, z2, s1, s2;
        slong prec1, prec2;

        acb_init(z1);
        acb_init(z2);
        acb_init(s1);
        acb_init(s2);

        acb_randtest(s1, state, 2 + n_randint(state, 100), 3);
        acb_randtest(z1, state, 2 + n_randint(state, 100), 3);
        acb_randtest(z2, state, 2 + n_randint(state, 100), 3);

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);

        acb_add(s2, s1, z1, prec2);
        acb_sub(s2, s2, z1, prec2);

        acb_dirichlet_eta(z1, s1, prec1);
        acb_dirichlet_eta(z2, s2, prec2);

        if (!acb_overlaps(z1, z2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("s1 = "); acb_printn(s1, 50, 0); flint_printf("\n\n");
            flint_printf("s2 = "); acb_printn(s2, 50, 0); flint_printf("\n\n");
            flint_printf("z1 = "); acb_printn(z1, 50, 0); flint_printf("\n\n");
            flint_printf("z2 = "); acb_printn(z2, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(s1);
        acb_clear(s2);
        acb_clear(z1);
        acb_clear(z2);
    }

    TEST_FUNCTION_END(state);
}
