/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"

TEST_FUNCTION_START(acb_vec_unit_roots, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong len, prec, order, k;
        acb_ptr vec;
        fmpq_t q;
        acb_t t;

        len = n_randint(state, 100);
        prec = 10 + n_randint(state, 200);

        switch (n_randint(state, 5))
        {
            case 0:
                order = 1 + len * n_randint(state, 6);
                break;
            case 1:
                order = 1 + n_randint(state, 4 * len);
                break;
            default:
                order = 1 + n_randint(state, len);
                break;
        }

        if (n_randint(state, 2))
            order = -order;

        vec = _acb_vec_init(len);
        acb_init(t);
        fmpq_init(q);

        _acb_vec_unit_roots(vec, order, len, prec);

        for (k = 0; k < len; k++)
        {
            if (order < 0)
                fmpq_set_si(q, -2 * k, -order);
            else
                fmpq_set_si(q, 2 * k, order);

            arb_sin_cos_pi_fmpq(acb_imagref(t), acb_realref(t), q, prec);

            if (!acb_overlaps(vec + k, t))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("n = %wu, order = %wd, k = %wd, prec = %wd\n\n", len, order, k, prec);
                flint_printf("vec = "); acb_printn(vec + k, 30, 0); flint_printf("\n\n");
                flint_printf("t = "); acb_printn(t, 30, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        _acb_vec_clear(vec, len);
        acb_clear(t);
        fmpq_clear(q);
    }

    TEST_FUNCTION_END(state);
}
