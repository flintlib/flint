/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_set_trunc_round, state)
{
    int iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_poly_t a, b, c, d, e;
        slong n, prec;

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);
        arb_poly_init(e);

        arb_poly_randtest(a, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);
        arb_poly_randtest(b, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);
        arb_poly_randtest(c, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);
        arb_poly_randtest(d, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);
        arb_poly_randtest(e, state, n_randint(state, 10), 2 + n_randint(state, 200), 10);

        n = n_randint(state, 10);
        prec = 2 + n_randint(state, 200);

        arb_poly_set_trunc(b, a, n);
        arb_poly_set_round(b, b, prec);

        arb_poly_set_round(c, a, prec);
        arb_poly_set_trunc(c, c, n);

        arb_poly_set_trunc_round(d, a, n, prec);

        arb_poly_set(e, a);
        arb_poly_set_trunc_round(e, e, n, prec);

        if (!arb_poly_equal(b, c) || !arb_poly_equal(c, d) || !arb_poly_equal(d, e))
        {
            flint_printf("FAIL\n\n");
            arb_poly_printd(a, 50), flint_printf("\n\n");
            arb_poly_printd(b, 50), flint_printf("\n\n");
            arb_poly_printd(c, 50), flint_printf("\n\n");
            arb_poly_printd(d, 50), flint_printf("\n\n");
            arb_poly_printd(e, 50), flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
        arb_poly_clear(e);
    }

    TEST_FUNCTION_END(state);
}
