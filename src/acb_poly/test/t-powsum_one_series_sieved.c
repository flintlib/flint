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

TEST_FUNCTION_START(acb_poly_powsum_one_series_sieved, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, a, q;
        acb_ptr z1, z2;
        slong i, n, len, prec;

        acb_init(s);
        acb_init(a);
        acb_init(q);

        if (n_randint(state, 2))
        {
            acb_randtest(s, state, 1 + n_randint(state, 200), 3);
        }
        else
        {
            arb_set_ui(acb_realref(s), 1);
            arb_mul_2exp_si(acb_realref(s), acb_realref(s), -1);
            arb_randtest(acb_imagref(s), state, 1 + n_randint(state, 200), 4);
        }

        acb_one(a);
        acb_one(q);

        prec = 2 + n_randint(state, 200);
        n = n_randint(state, 100);
        len = 1 + n_randint(state, 10);

        z1 = _acb_vec_init(len);
        z2 = _acb_vec_init(len);

        _acb_poly_powsum_series_naive(z1, s, a, q, n, len, prec);
        _acb_poly_powsum_one_series_sieved(z2, s, n, len, prec);

        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(z1 + i, z2 + i))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("iter = %wd\n", iter);
                flint_printf("n = %wd, prec = %wd, len = %wd, i = %wd\n\n", n, prec, len, i);
                flint_printf("s = "); acb_printd(s, prec / 3.33); flint_printf("\n\n");
                flint_printf("a = "); acb_printd(a, prec / 3.33); flint_printf("\n\n");
                flint_printf("z1 = "); acb_printd(z1 + i, prec / 3.33); flint_printf("\n\n");
                flint_printf("z2 = "); acb_printd(z2 + i, prec / 3.33); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(s);
        acb_clear(q);
        _acb_vec_clear(z1, len);
        _acb_vec_clear(z2, len);
    }

    TEST_FUNCTION_END(state);
}
