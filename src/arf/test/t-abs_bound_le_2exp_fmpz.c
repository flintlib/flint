/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_abs_bound_le_2exp_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y;
        fmpz_t b;
        int cmp1, cmp2;

        arf_init(x);
        arf_init(y);
        fmpz_init(b);

        arf_randtest_not_zero(x, state, 2 + n_randint(state, 1000), 100);
        arf_abs_bound_le_2exp_fmpz(b, x);

        arf_one(y);
        arf_mul_2exp_fmpz(y, y, b);

        cmp1 = (arf_cmpabs(x, y) <= 0);

        arf_mul_2exp_si(y, y, -1);

        cmp2 = (arf_cmpabs(y, x) < 0);

        arf_mul_2exp_si(y, y, 1);

        if (!cmp1 || !cmp2)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
            flint_printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        fmpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
