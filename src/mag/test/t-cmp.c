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

TEST_FUNCTION_START(mag_cmp, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y;
        mag_t xb, yb;
        int c1, c2;

        arf_init(x);
        arf_init(y);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 100);
        mag_randtest_special(yb, state, 100);

        arf_set_mag(x, xb);
        arf_set_mag(y, yb);

        c1 = arf_cmp(x, y);
        c2 = mag_cmp(xb, yb);

        if (c1 != c2)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("xb = "); mag_print(xb); flint_printf("\n\n");
            flint_printf("yb = "); mag_print(yb); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);

        mag_clear(xb);
        mag_clear(yb);
    }

    TEST_FUNCTION_END(state);
}
