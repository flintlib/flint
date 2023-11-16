/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_contains_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x;
        fmpq_t y;
        fmpz_t a, b, t;
        int c1, c2;
        slong shift;

        arb_init(x);
        fmpq_init(y);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(t);

        arb_randtest_special(x, state, 1 + n_randint(state, 500), 12);

        if (n_randint(state, 2) && arb_is_finite(x))
        {
            arb_get_rand_fmpq(y, state, x, 1 + n_randint(state, 500));

            fmpz_add_si(fmpq_numref(y), fmpq_numref(y), n_randint(state, 3) - 1);
        }
        else
            fmpq_randtest(y, state, 1 + n_randint(state, 500));

        c1 = arb_contains_fmpq(x, y);

        if (arb_is_finite(x))
        {
            arb_get_interval_fmpz_2exp(a, b, t, x);
            shift = fmpz_get_si(t);

            fmpz_mul(a, a, fmpq_denref(y));
            fmpz_mul(b, b, fmpq_denref(y));

            fmpz_set(t, fmpq_numref(y));

            if (shift >= 0)
            {
                fmpz_mul_2exp(a, a, shift);
                fmpz_mul_2exp(b, b, shift);
            }
            else
            {
                fmpz_mul_2exp(t, t, -shift);
            }

            c2 = (fmpz_cmp(a, t) <= 0) && (fmpz_cmp(t, b) <= 0);
        }
        else
        {
            fmpz_tdiv_q(t, fmpq_numref(y), fmpq_denref(y));
            c2 = arb_contains_fmpz(x, t);
        }

        if (c1 != c2)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("c1 = %d, c2 = %d\n\n", c1, c2);
            flint_abort();
        }

        arb_clear(x);
        fmpq_clear(y);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
