/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"
#include "mag.h"
#include "long_extras.h"

TEST_FUNCTION_START(mag_cmp_2exp_si, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x;
        mag_t xb;
        slong y;
        int c1, c2;

        arf_init(x);
        mag_init(xb);

        mag_randtest_special(xb, state, 100);
        y = z_randtest(state);

        arf_set_mag(x, xb);

        c1 = arf_cmp_2exp_si(x, y);
        c2 = mag_cmp_2exp_si(xb, y);

        if (c1 != c2)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = %wd", y);  flint_printf("\n\n");
            flint_printf("xb = "); mag_print(xb); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        mag_clear(xb);
    }

    TEST_FUNCTION_END(state);
}
