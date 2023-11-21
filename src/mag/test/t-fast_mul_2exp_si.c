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

TEST_FUNCTION_START(mag_fast_mul_2exp_si, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z;
        mag_t xb, yb;
        slong e;

        arf_init(x);
        arf_init(y);
        arf_init(z);

        mag_init(xb);
        mag_init(yb);

        mag_randtest(xb, state, 15);
        e = n_randint(state, 10000) - n_randint(state, 10000);
        arf_set_mag(x, xb);

        mag_fast_mul_2exp_si(yb, xb, e);

        arf_mul_2exp_si(y, x, e);

        arf_set_mag(z, yb);

        MAG_CHECK_BITS(yb)

        if (!arf_equal(z, y))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); arf_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); arf_printd(z, 15); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);

        mag_clear(xb);
        mag_clear(yb);
    }

    TEST_FUNCTION_END(state);
}
