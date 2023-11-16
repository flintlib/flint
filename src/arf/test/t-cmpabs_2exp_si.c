/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_cmpabs_2exp_si, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong bits, e;
        arf_t x, y;
        int cmp1, cmp2;

        bits = 2 + n_randint(state, 1000);
        e = n_randtest(state);

        arf_init(x);
        arf_init(y);

        if (iter % 10 == 0)
            arf_set_ui_2exp_si(x, 1, e);
        else
            arf_randtest_special(x, state, bits, 100);

        arf_set_ui_2exp_si(y, 1, e);

        cmp1 = arf_cmpabs(x, y);
        cmp2 = arf_cmpabs_2exp_si(x, e);

        if (cmp1 != cmp2)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("cmp1 = %d, cmp2 = %d\n\n", cmp1, cmp2);
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
    }

    TEST_FUNCTION_END(state);
}
