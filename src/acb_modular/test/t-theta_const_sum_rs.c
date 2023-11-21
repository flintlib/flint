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
#include "acb_modular.h"

TEST_FUNCTION_START(acb_modular_theta_const_sum_rs, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t q, t2a, t2b, t3a, t3b, t4a, t4b;
        slong prec1, prec2, N;

        acb_init(q);
        acb_init(t2a);
        acb_init(t2b);
        acb_init(t3a);
        acb_init(t3b);
        acb_init(t4a);
        acb_init(t4b);

        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);
        N = n_randint(state, 300);

        acb_randtest(q, state, prec1, 3);

        acb_randtest(t2a, state, prec1, 3);
        acb_randtest(t2b, state, prec1, 3);
        acb_randtest(t3a, state, prec1, 3);
        acb_randtest(t3b, state, prec1, 3);
        acb_randtest(t4a, state, prec1, 3);
        acb_randtest(t4b, state, prec1, 3);

        acb_modular_theta_const_sum_basecase(t2a, t3a, t4a, q, N, prec1);
        acb_modular_theta_const_sum_rs(t2b, t3b, t4b, q, N, prec2);

        if (!acb_overlaps(t2a, t2b) || !acb_overlaps(t3a, t3b) || !acb_overlaps(t4a, t4b))
        {
            flint_printf("FAIL (overlap)  iter = %wd\n", iter);
            flint_printf("N = %wd\n", N);
            flint_printf("q = "); acb_printd(q, 50); flint_printf("\n\n");
            flint_printf("t2a = "); acb_printd(t2a, 50); flint_printf("\n\n");
            flint_printf("t2b = "); acb_printd(t2b, 50); flint_printf("\n\n");
            flint_printf("t3a = "); acb_printd(t3a, 50); flint_printf("\n\n");
            flint_printf("t3b = "); acb_printd(t3b, 50); flint_printf("\n\n");
            flint_printf("t4a = "); acb_printd(t4a, 50); flint_printf("\n\n");
            flint_printf("t4b = "); acb_printd(t4b, 50); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(q);
        acb_clear(t2a);
        acb_clear(t2b);
        acb_clear(t3a);
        acb_clear(t3b);
        acb_clear(t4a);
        acb_clear(t4b);
    }

    TEST_FUNCTION_END(state);
}
