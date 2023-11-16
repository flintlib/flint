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

TEST_FUNCTION_START(acb_poly_exp_pi_i_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, prec;
        acb_poly_t a, b, c;
        acb_t t;

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_init(t);

        prec = 2 + n_randint(state, 100);

        acb_poly_randtest(a, state, 1 + n_randint(state, 20), prec, 5);
        acb_poly_randtest(b, state, 1 + n_randint(state, 20), prec, 5);
        acb_poly_randtest(c, state, 1 + n_randint(state, 20), prec, 5);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        if (n_randint(state, 2) == 0)
            acb_poly_exp_pi_i_series(b, a, m, prec);
        else
        {
            acb_poly_set(b, a);
            acb_poly_exp_pi_i_series(b, b, m, prec);
        }

        acb_const_pi(t, prec);
        acb_mul_onei(t, t);
        acb_poly_scalar_mul(c, a, t, prec);
        acb_poly_exp_series(c, c, n, prec);

        acb_poly_truncate(b, FLINT_MIN(m, n));
        acb_poly_truncate(c, FLINT_MIN(m, n));

        if (!acb_poly_overlaps(b, c))
        {
            flint_printf("FAIL\n\n");
            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
