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

TEST_FUNCTION_START(arb_sin_pi_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t s1, s2;
        fmpq_t x;
        slong prec;

        prec = 2 + n_randint(state, 5000);

        arb_init(s1);
        arb_init(s2);
        fmpq_init(x);

        fmpq_randtest(x, state, 1 + n_randint(state, 200));

        arb_sin_pi_fmpq(s1, x, prec);

        arb_const_pi(s2, prec);
        arb_mul_fmpz(s2, s2, fmpq_numref(x), prec);
        arb_div_fmpz(s2, s2, fmpq_denref(x), prec);
        arb_sin(s2, s2, prec);

        if (!arb_overlaps(s1, s2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printd(s1, 15); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printd(s2, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(s1);
        arb_clear(s2);
        fmpq_clear(x);
    }

    TEST_FUNCTION_END(state);
}
