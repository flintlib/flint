/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_pfq, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_ptr a, b;
        acb_t z, t, res1, res2;
        slong i, p, q, extend, prec1, prec2;
        int regularized1, regularized2;

        p = n_randint(state, 3);
        q = n_randint(state, 3);
        extend = n_randint(state, 3);

        a = _acb_vec_init(p + extend);
        b = _acb_vec_init(q + extend);
        acb_init(z);
        acb_init(t);
        acb_init(res1);
        acb_init(res2);

        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);
        regularized1 = n_randint(state, 2);
        regularized2 = n_randint(state, 2);

        for (i = 0; i < p; i++)
        {
            if (n_randint(state, 2))
                acb_one(a + i);
            else
                acb_randtest_param(a + i, state, 1 + n_randint(state, 200), 3);
        }

        for (i = 0; i < q; i++)
            acb_randtest_param(b + i, state, 1 + n_randint(state, 200), 3);

        for (i = 0; i < extend; i++)
        {
            acb_randtest_param(a + p + i, state, 1 + n_randint(state, 200), 3);
            acb_set(b + q + i, a + p + i);
        }

        acb_randtest_param(z, state, 1 + n_randint(state, 200), 3);

        acb_hypgeom_pfq(res1, a, p, b, q, z, regularized1, prec1);
        acb_hypgeom_pfq(res2, a, p + extend, b, q + extend, z, regularized2, prec2);

        if (regularized1 && !regularized2)
        {
            for (i = 0; i < q; i++)
            {
                acb_rgamma(t, b + i, prec2);
                acb_mul(res2, res2, t, prec2);
            }
        }
        else if (regularized2)
        {
            if (!regularized1)
            {
                for (i = 0; i < q; i++)
                {
                    acb_rgamma(t, b + i, prec2);
                    acb_mul(res1, res1, t, prec2);
                }
            }

            for (i = 0; i < extend; i++)
            {
                acb_rgamma(t, b + q + i, prec2);
                acb_mul(res1, res1, t, prec2);
            }
        }

        if (!acb_overlaps(res1, res2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("p = %wd, q = %wd, regularized1 = %d, regularized2 = %d\n\n",
                p, q, regularized1, regularized2);

            for (i = 0; i < p + extend; i++)
            {
                flint_printf("a[%wd] = ", i);
                acb_printd(a + i, 30);
                flint_printf("\n\n");
            }

            for (i = 0; i < q + extend; i++)
            {
                flint_printf("b[%wd] = ", i);
                acb_printd(b + i, 30);
                flint_printf("\n\n");
            }

            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("res1 = "); acb_printd(res1, 30); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 30); flint_printf("\n\n");
            flint_abort();
        }

        _acb_vec_clear(a, p + extend);
        _acb_vec_clear(b, q + extend);
        acb_clear(z);
        acb_clear(t);
        acb_clear(res1);
        acb_clear(res2);
    }

    TEST_FUNCTION_END(state);
}
