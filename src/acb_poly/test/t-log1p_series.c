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

TEST_FUNCTION_START(acb_poly_log1p_series, state)
{
    slong iter;

    /* compare with log_series */
    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t a, b, c;
        slong m, n, prec;

        prec = 2 + n_randint(state, 200);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);

        acb_poly_randtest(a, state, 1 + n_randint(state, 10), prec, 5);
        acb_poly_randtest(b, state, 1 + n_randint(state, 10), prec, 5);
        acb_poly_randtest(c, state, 1 + n_randint(state, 10), prec, 5);

        if (n_randint(state, 2))
            acb_poly_log1p_series(b, a, n, prec);
        else
        {
            acb_poly_set(b, a);
            acb_poly_log1p_series(b, b, n, prec);
        }

        acb_poly_add_si(c, a, 1, prec);
        acb_poly_log_series(c, c, m, prec);

        acb_poly_truncate(b, FLINT_MIN(m, n));
        acb_poly_truncate(c, FLINT_MIN(m, n));

        if (!acb_poly_overlaps(b, c))
        {
            flint_printf("FAIL\n\n");
            flint_printf("m = %wd, n = %wd\n\n", m, n);
            flint_printf("a = "); acb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_poly_printd(c, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(a);
        acb_poly_clear(b);
        acb_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
