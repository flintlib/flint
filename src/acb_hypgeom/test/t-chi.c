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
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_chi, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z0, z1, w0, w1;
        slong prec0, prec1;

        acb_init(z0);
        acb_init(z1);
        acb_init(w0);
        acb_init(w1);

        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);

        acb_randtest(z0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));

        acb_set(z1, z0);
        if (n_randint(state, 2))
        {
            acb_add(z1, z1, w0, prec0);
            acb_sub(z1, z1, w0, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_chi_2f3(w0, z0, prec0);
                break;
            case 1:
                acb_hypgeom_chi_asymp(w0, z0, prec0);
                break;
            default:
                acb_hypgeom_chi(w0, z0, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_chi_2f3(w1, z1, prec1);
                break;
            case 1:
                acb_hypgeom_chi_asymp(w1, z1, prec1);
                break;
            default:
                acb_hypgeom_chi(w1, z1, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("z0 = "); acb_printd(z0, 30); flint_printf("\n\n");
            flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z0);
        acb_clear(z1);
        acb_clear(w0);
        acb_clear(w1);
    }

    TEST_FUNCTION_END(state);
}
