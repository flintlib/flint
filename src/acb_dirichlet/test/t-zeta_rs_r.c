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

TEST_FUNCTION_START(acb_dirichlet_zeta_rs_r, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z1, z2, s1, s2;
        slong prec1, prec2, K1, K2;

        acb_init(z1);
        acb_init(z2);
        acb_init(s1);
        acb_init(s2);

        if (n_randint(state, 2))
            arb_set_d(acb_realref(s1), 0.5);
        else
            arb_randtest(acb_realref(s1), state, 2 + n_randint(state, 100), 2);

        arb_randtest(acb_imagref(s1), state, 2 + n_randint(state, 100), 2);
        arb_abs(acb_imagref(s1), acb_imagref(s1));
        arb_add_ui(acb_imagref(s1), acb_imagref(s1), n_randtest(state) % 1000000, 100);

        acb_randtest(z1, state, 2 + n_randint(state, 100), 3);
        acb_randtest(z2, state, 2 + n_randint(state, 100), 3);

        prec1 = 2 + n_randint(state, 150);
        prec2 = 2 + n_randint(state, 150);

        K1 = n_randint(state, 10);
        K2 = n_randint(state, 10);

        if (n_randint(state, 2))
        {
            mag_zero(arb_radref(acb_realref(s1)));
            mag_zero(arb_radref(acb_imagref(s1)));
        }

        if (n_randint(state, 2))
        {
            acb_set(s2, s1);
        }
        else
        {
            acb_add(s2, s1, z1, prec2);
            acb_sub(s2, s2, z1, prec2);
        }

        acb_dirichlet_zeta_rs_r(z1, s1, K1, prec1);
        acb_dirichlet_zeta_rs_r(z2, s2, K2, prec2);

        if (!acb_overlaps(z1, z2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd\n", iter);
            flint_printf("K1 = %wd, K2 = %wd\n\n", K1, K2);
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
