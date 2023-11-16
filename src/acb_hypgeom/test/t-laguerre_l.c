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

TEST_FUNCTION_START(acb_hypgeom_laguerre_l, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t n, m, n1, m1, z, res1, res2, res3, s;
        slong prec;

        acb_init(n);
        acb_init(m);
        acb_init(n1);
        acb_init(m1);
        acb_init(z);
        acb_init(res1);
        acb_init(res2);
        acb_init(res3);
        acb_init(s);

        prec = 2 + n_randint(state, 200);

        if (n_randint(state, 2))
        {
            acb_set_si(n, n_randint(state, 20) - 10);
            acb_set_si(m, n_randint(state, 20) - 10);
        }
        else
        {
            acb_randtest_param(n, state, 1 + n_randint(state, 400), 10);
            acb_randtest_param(m, state, 1 + n_randint(state, 400), 10);
        }

        acb_randtest_param(z, state, 1 + n_randint(state, 400), 10);

        acb_sub_ui(n1, n, 1, prec);
        acb_sub_ui(m1, m, 1, prec);

        acb_hypgeom_laguerre_l(res1, n, m, z, prec);
        acb_hypgeom_laguerre_l(res2, n1, m, z, 2 + n_randint(state, 200));
        acb_hypgeom_laguerre_l(res3, n, m1, z, 2 + n_randint(state, 200));

        acb_add(s, res2, res3, prec);

        if (acb_is_finite(res1) && acb_is_finite(s) && !acb_overlaps(res1, s))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("iter = %wd\n\n", iter);
            flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
            flint_printf("m = "); acb_printd(m, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("res1 = "); acb_printd(res1, 30); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 30); flint_printf("\n\n");
            flint_printf("res3 = "); acb_printd(res3, 30); flint_printf("\n\n");
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(n);
        acb_clear(m);
        acb_clear(n1);
        acb_clear(m1);
        acb_clear(z);
        acb_clear(res1);
        acb_clear(res2);
        acb_clear(res3);
        acb_clear(s);
    }

    TEST_FUNCTION_END(state);
}
