/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_zeta_cpx_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, a;
        acb_ptr z1, z2;
        slong i, len, prec1, prec2;
        int deflate;

        acb_init(s);
        acb_init(a);

        if (n_randint(state, 2))
        {
            acb_randtest(s, state, 1 + n_randint(state, 300), 3);
        }
        else
        {
            arb_set_ui(acb_realref(s), 1);
            arb_mul_2exp_si(acb_realref(s), acb_realref(s), -1);
            arb_randtest(acb_imagref(s), state, 1 + n_randint(state, 300), 4);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_randtest(a, state, 1 + n_randint(state, 300), 3);
                break;
            case 1:
                arb_randtest(acb_realref(a), state, 1 + n_randint(state, 300), 3);
                break;
            case 2:
                acb_one(a);
                break;
        }

        prec1 = 2 + n_randint(state, 300);
        prec2 = prec1 + 30;
        len = 1 + n_randint(state, 20);

        deflate = n_randint(state, 2);

        z1 = _acb_vec_init(len);
        z2 = _acb_vec_init(len);

        _acb_poly_zeta_cpx_series(z1, s, a, deflate, len, prec1);
        _acb_poly_zeta_cpx_series(z2, s, a, deflate, len, prec2);

        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(z1 + i, z2 + i))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("iter = %wd\n", iter);
                flint_printf("deflate = %d, len = %wd, i = %wd\n\n", deflate, len, i);
                flint_printf("s = "); acb_printd(s, prec1 / 3.33); flint_printf("\n\n");
                flint_printf("a = "); acb_printd(a, prec1 / 3.33); flint_printf("\n\n");
                flint_printf("z1 = "); acb_printd(z1 + i, prec1 / 3.33); flint_printf("\n\n");
                flint_printf("z2 = "); acb_printd(z2 + i, prec2 / 3.33); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(s);
        _acb_vec_clear(z1, len);
        _acb_vec_clear(z2, len);
    }

    TEST_FUNCTION_END(state);
}
