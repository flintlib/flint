/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_floor, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t x;
        fmpz_t y;
        slong prec;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpz_init(y);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);
        prec = 2 + n_randint(state, 200);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));

        arb_floor(b, a, prec);
        fmpz_fdiv_q(y, fmpq_numref(x), fmpq_denref(x));

        if (!arb_contains_fmpz(b, y))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arb_floor(a, a, prec);

        if (!arb_equal(a, b))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpz_clear(y);
    }

    TEST_FUNCTION_END(state);
}
