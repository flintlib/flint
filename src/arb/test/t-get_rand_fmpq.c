/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_get_rand_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x;
        fmpq_t q;

        arb_init(x);
        fmpq_init(q);

        arb_randtest(x, state, 1 + n_randint(state, 200), 10);
        arb_get_rand_fmpq(q, state, x, 1 + n_randint(state, 200));

        if (!arb_contains_fmpq(x, q) || !fmpq_is_canonical(q))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("q = "); fmpq_print(q); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        fmpq_clear(q);
    }

    TEST_FUNCTION_END(state);
}
