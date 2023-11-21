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

TEST_FUNCTION_START(acb_hypgeom_hermite_h, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t n, n1, n2, z, res1, res2, res3, s;
        slong prec1, prec2, prec3;

        acb_init(n);
        acb_init(n1);
        acb_init(n2);
        acb_init(z);
        acb_init(res1);
        acb_init(res2);
        acb_init(res3);
        acb_init(s);

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);
        prec3 = 2 + n_randint(state, 300);

        acb_randtest_param(n, state, 1 + n_randint(state, 400), 10);
        acb_sub_ui(n1, n, 1, prec1);
        acb_sub_ui(n2, n, 2, prec1);
        acb_randtest_param(z, state, 1 + n_randint(state, 400), 10);
        acb_randtest_param(res1, state, 1 + n_randint(state, 400), 10);

        acb_hypgeom_hermite_h(res1, n, z, prec1);
        acb_hypgeom_hermite_h(res2, n1, z, prec2);
        acb_hypgeom_hermite_h(res3, n2, z, prec3);

        acb_mul(s, res2, z, prec1);
        acb_submul(s, res3, n1, prec1);
        acb_mul_2exp_si(s, s, 1);

        if (!acb_overlaps(res1, s))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("iter = %wd\n\n", iter);
            flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("res1 = "); acb_printd(res1, 30); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 30); flint_printf("\n\n");
            flint_printf("res3 = "); acb_printd(res3, 30); flint_printf("\n\n");
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(n);
        acb_clear(n1);
        acb_clear(n2);
        acb_clear(z);
        acb_clear(res1);
        acb_clear(res2);
        acb_clear(res3);
        acb_clear(s);
    }

    TEST_FUNCTION_END(state);
}
