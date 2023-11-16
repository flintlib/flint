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

TEST_FUNCTION_START(acb_poly_exp_series_basecase, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, bits1, bits2, bits3;
        acb_poly_t a, b, c, d;

        bits1 = 2 + n_randint(state, 100);
        bits2 = 2 + n_randint(state, 100);
        bits3 = 2 + n_randint(state, 100);

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        acb_poly_randtest(a, state, m, bits1, 5);
        acb_poly_randtest(b, state, m, bits1, 5);

        /* check exp(a+b) = exp(a) exp(b) */
        acb_poly_exp_series_basecase(c, a, n, bits2);
        acb_poly_exp_series_basecase(d, b, n, bits2);
        acb_poly_mullow(c, c, d, n, bits2);

        acb_poly_add(d, a, b, bits3);
        acb_poly_exp_series_basecase(d, d, n, bits3);   /* also aliasing test */

        if (!acb_poly_overlaps(c, d))
        {
            flint_printf("FAIL\n\n");

            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); acb_poly_printd(d, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
