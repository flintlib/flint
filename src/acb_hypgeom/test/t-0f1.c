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

TEST_FUNCTION_START(acb_hypgeom_0f1, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, z, w1, w2;
        slong prec1, prec2;
        int regularized, ebits;

        acb_init(a);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        prec1 = 2 + n_randint(state, 700);
        prec2 = 2 + n_randint(state, 700);

        if (n_randint(state, 5) == 0)
            ebits = 100;
        else
            ebits = 10;

        acb_randtest_param(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest_param(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        regularized = n_randint(state, 2);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_0f1_asymp(w1, a, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_0f1_direct(w1, a, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_0f1(w1, a, z, regularized, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_0f1_asymp(w2, a, z, regularized, prec2);
                break;
            case 1:
                acb_hypgeom_0f1_direct(w2, a, z, regularized, prec2);
                break;
            default:
                acb_hypgeom_0f1(w2, a, z, regularized, prec2);
        }

        if (!acb_overlaps(w1, w2))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
    }

    TEST_FUNCTION_END(state);
}
