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

TEST_FUNCTION_START(acb_hypgeom_chebyshev_u, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t n, z, a, b, c, t, res1, res2;
        slong prec;

        acb_init(n);
        acb_init(z);
        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(t);
        acb_init(res1);
        acb_init(res2);

        prec = 2 + n_randint(state, 300);

        acb_randtest_param(n, state, 1 + n_randint(state, 400), 1 + n_randint(state, 10));
        acb_randtest_param(z, state, 1 + n_randint(state, 400), 1 + n_randint(state, 10));
        acb_randtest_param(res1, state, 1 + n_randint(state, 400), 1 + n_randint(state, 10));
        acb_randtest_param(res2, state, 1 + n_randint(state, 400), 1 + n_randint(state, 10));

        acb_hypgeom_chebyshev_u(res1, n, z, prec);

        acb_neg(a, n);
        acb_add_ui(b, n, 2, prec);
        acb_set_ui(c, 3);
        acb_mul_2exp_si(c, c, -1);
        acb_sub_ui(t, z, 1, prec);
        acb_neg(t, t);
        acb_mul_2exp_si(t, t, -1);
        acb_hypgeom_2f1(res2, a, b, c, t, 0, 2 + n_randint(state, 300));
        acb_add_ui(t, n, 1, prec);
        acb_mul(res2, res2, t, prec);

        if (!acb_overlaps(res1, res2))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("iter = %wd\n\n", iter);
            flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("res1 = "); acb_printd(res1, 30); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(n);
        acb_clear(z);
        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(t);
        acb_clear(res1);
        acb_clear(res2);
    }

    TEST_FUNCTION_END(state);
}
