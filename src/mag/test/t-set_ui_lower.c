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

TEST_FUNCTION_START(mag_set_ui_lower, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t a, b, c;
        mag_t m;
        ulong x;

        arf_init(a);
        arf_init(b);
        arf_init(c);
        mag_init(m);

        x = n_randtest(state);

        arf_set_ui(a, x);
        mag_set_ui_lower(m, x);

        arf_set_mag(b, m);

        arf_set(c, a);
        arf_mul_ui(c, c, 1023, MAG_BITS, ARF_RND_DOWN);
        arf_mul_2exp_si(c, c, -10);

        MAG_CHECK_BITS(m)

        if (!(arf_cmpabs(c, b) <= 0 && arf_cmpabs(b, a) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = %wu\n\n", x);
            flint_printf("a = "); arf_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("c = "); arf_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(a);
        arf_clear(b);
        arf_clear(c);
        mag_clear(m);
    }

    TEST_FUNCTION_END(state);
}
