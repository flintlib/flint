/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_evaluate_rectangular, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_poly_t f;
        arb_t x, y1, y2;

        arb_init(x);
        arb_init(y1);
        arb_init(y2);
        arb_poly_init(f);

        arb_randtest(x, state, 2 + n_randint(state, 1000), 5);
        arb_poly_randtest(f, state, 2 + n_randint(state, 100), 2 + n_randint(state, 1000), 5);

        arb_poly_evaluate_rectangular(y1, f, x, 2 + n_randint(state, 1000));
        arb_poly_evaluate_horner(y2, f, x, 2 + n_randint(state, 1000));

        if (!arb_overlaps(y1, y2))
        {
            flint_printf("FAIL\n\n");
            flint_printf("f = "); arb_poly_printd(f, 15); flint_printf("\n\n");
            flint_printf("x = "); arb_printd(x, 15); flint_printf("\n\n");
            flint_printf("y1 = "); arb_printd(y1, 15); flint_printf("\n\n");
            flint_printf("y2 = "); arb_printd(y2, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_clear(f);
        arb_clear(x);
        arb_clear(y1);
        arb_clear(y2);
    }

    TEST_FUNCTION_END(state);
}
