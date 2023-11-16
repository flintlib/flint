/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_set_trunc_round, state)
{
    int iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t a, b, c, d, e;
        slong n, prec;

        acb_poly_init(a);
        acb_poly_init(b);
        acb_poly_init(c);
        acb_poly_init(d);
        acb_poly_init(e);

        acb_poly_randtest(a, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);
        acb_poly_randtest(b, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);
        acb_poly_randtest(c, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);
        acb_poly_randtest(d, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);
        acb_poly_randtest(e, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);

        n = n_randint(state, 10);
        prec = 2 + n_randint(state, 200);

        acb_poly_set_trunc(b, a, n);
        acb_poly_set_round(b, b, prec);

        acb_poly_set_round(c, a, prec);
        acb_poly_set_trunc(c, c, n);

        acb_poly_set_trunc_round(d, a, n, prec);

        acb_poly_set(e, a);
        acb_poly_set_trunc_round(e, e, n, prec);

        if (!acb_poly_equal(b, c) || !acb_poly_equal(c, d) || !acb_poly_equal(d, e))
        {
            flint_printf("FAIL\n\n");
            acb_poly_printd(a, 50), flint_printf("\n\n");
            acb_poly_printd(b, 50), flint_printf("\n\n");
            acb_poly_printd(c, 50), flint_printf("\n\n");
            acb_poly_printd(d, 50), flint_printf("\n\n");
            acb_poly_printd(e, 50), flint_printf("\n\n");
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
