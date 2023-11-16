/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_cos_pi_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t c1, c2;
        fmpq_t x;
        slong prec;

        prec = 2 + n_randint(state, 5000);

        arb_init(c1);
        arb_init(c2);
        fmpq_init(x);

        fmpq_randtest(x, state, 1 + n_randint(state, 200));

        arb_cos_pi_fmpq(c1, x, prec);

        arb_const_pi(c2, prec);
        arb_mul_fmpz(c2, c2, fmpq_numref(x), prec);
        arb_div_fmpz(c2, c2, fmpq_denref(x), prec);
        arb_cos(c2, c2, prec);

        if (!arb_overlaps(c1, c2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("c1 = "); arb_printd(c1, 15); flint_printf("\n\n");
            flint_printf("c2 = "); arb_printd(c2, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(c1);
        arb_clear(c2);
        fmpq_clear(x);
    }

    TEST_FUNCTION_END(state);
}
