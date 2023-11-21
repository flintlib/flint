/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_is_int_2exp_si, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y;
        fmpz_t t;
        slong e;
        int res1, res2;

        arf_init(x);
        arf_init(y);
        fmpz_init(t);

        arf_randtest_special(x, state, 2000, 100);
        e = n_randtest(state);
        arf_mul_2exp_si(y, x, e);

        res1 = arf_is_int(x);
        res2 = arf_is_int_2exp_si(y, e);

        if (res1 != res2)
        {
            flint_printf("FAIL! (1)\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("res1 = %d, res2 = %d\n\n", res1, res2);
            flint_abort();
        }

        if (res1)
        {
            if (n_randint(state, 2))
                arf_floor(y, x);
            else
                arf_ceil(y, x);

            if (!arf_equal(x, y) || !arf_is_finite(x))
            {
                flint_printf("FAIL! (2)\n");
                flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                flint_printf("res1 = %d\n\n", res1);
                flint_abort();
            }
        }

        if (arf_is_finite(x) && !arf_is_zero(x))
        {
            arf_bot(t, x);
            fmpz_neg(t, t);
            arf_mul_2exp_fmpz(x, x, t);
            res1 = arf_is_int(x);
            arf_mul_2exp_si(y, x, -1);
            res2 = arf_is_int(y);

            if (!arf_is_int(x) || arf_is_int(y))
            {
                flint_printf("FAIL! (3)\n");
                flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                flint_printf("res1 = %d, res2 = %d\n\n", res1, res2);
                flint_abort();
            }
        }

        arf_clear(x);
        arf_clear(y);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
