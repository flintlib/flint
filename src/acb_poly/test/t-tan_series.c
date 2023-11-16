/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_tan_series, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, rbits1, rbits2;
        acb_poly_t a, b, c, d, e;

        rbits1 = 2 + n_randint(state, 100);
        rbits2 = 2 + n_randint(state, 100);

        if (n_randint(state, 10) == 0)
        {
            m = 1 + n_randint(state, 50);
            n = 1 + n_randint(state, 50);
        }
        else
        {
            m = 1 + n_randint(state, 20);
            n = 1 + n_randint(state, 20);
        }

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);
        acb_poly_init(e);

        acb_poly_randtest(a, state, m, rbits1, 10);
        acb_poly_tan_series(b, a, n, rbits2);

        /* check tan(x) = 2*tan(x/2)/(1-tan(x/2)^2) */
        acb_poly_scalar_mul_2exp_si(c, a, -1);
        acb_poly_tan_series(c, c, n, rbits2);
        acb_poly_mullow(d, c, c, n, rbits2);
        acb_poly_one(e);
        acb_poly_sub(e, e, d, rbits2);
        acb_poly_div_series(c, c, e, n, rbits2);
        acb_poly_scalar_mul_2exp_si(c, c, 1);

        if (!acb_poly_overlaps(b, c))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits2 = %wd\n", rbits2);

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_tan_series(a, a, n, rbits2);
        if (!acb_poly_equal(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
        acb_poly_clear(e);
    }

    TEST_FUNCTION_END(state);
}
