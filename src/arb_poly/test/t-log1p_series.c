/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_log1p_series, state)
{
    slong iter;

    /* compare with log_series */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_poly_t a, b, c;
        slong m, n, prec;

        prec = 2 + n_randint(state, 200);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);

        arb_poly_randtest(a, state, 1 + n_randint(state, 10), prec, 5);
        arb_poly_randtest(b, state, 1 + n_randint(state, 10), prec, 5);
        arb_poly_randtest(c, state, 1 + n_randint(state, 10), prec, 5);

        if (n_randint(state, 2))
            arb_poly_log1p_series(b, a, n, prec);
        else
        {
            arb_poly_set(b, a);
            arb_poly_log1p_series(b, b, n, prec);
        }

        arb_poly_add_si(c, a, 1, prec);
        arb_poly_log_series(c, c, m, prec);

        arb_poly_truncate(b, FLINT_MIN(m, n));
        arb_poly_truncate(c, FLINT_MIN(m, n));

        if (!arb_poly_overlaps(b, c))
        {
            flint_printf("FAIL\n\n");
            flint_printf("m = %wd, n = %wd\n\n", m, n);
            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
    }

    TEST_FUNCTION_END(state);
}
