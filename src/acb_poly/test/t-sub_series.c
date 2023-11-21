/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_sub_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t a, b, c, d;
        slong len, prec;

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);

        acb_poly_randtest(a, state, 1 + n_randint(state, 10), 100, 10);
        acb_poly_randtest(b, state, 1 + n_randint(state, 10), 100, 10);
        acb_poly_randtest(c, state, 1 + n_randint(state, 10), 100, 10);
        acb_poly_randtest(d, state, 1 + n_randint(state, 10), 100, 10);
        prec = 2 + n_randint(state, 100);
        len = n_randint(state, 10);

        acb_poly_sub_series(c, a, b, len, prec);
        acb_poly_sub(d, a, b, prec);
        acb_poly_truncate(d, len);

        if (!acb_poly_equal(c, d))
        {
            flint_printf("FAIL\n\n");
            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_set(d, a);
        acb_poly_sub_series(d, d, b, len, prec);
        if (!acb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 1)\n\n");
            flint_abort();
        }

        acb_poly_set(d, b);
        acb_poly_sub_series(d, a, d, len, prec);
        if (!acb_poly_equal(d, c))
        {
            flint_printf("FAIL (aliasing 2)\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
        acb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
