/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_rising_ui_get_mag, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, y, z;
        mag_t b;
        ulong n;

        acb_init(x);
        acb_init(y);
        acb_init(z);
        mag_init(b);

        acb_randtest(x, state, 1 + n_randint(state, 400), 1 + n_randint(state, 100));
        n = n_randint(state, 80);

        acb_rising_ui(y, x, n, 2 + n_randint(state, 400));
        acb_rising_ui_get_mag(b, x, n);
        acb_zero(z);
        acb_add_error_mag(z, b);

        if (!acb_overlaps(y, z))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu\n", n);
            flint_printf("x = "); acb_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); acb_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(z);
        mag_clear(b);
    }

    TEST_FUNCTION_END(state);
}
