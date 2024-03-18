/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"
#include "mag.h"

TEST_FUNCTION_START(mag_bin_uiui, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y;
        fmpz_t f;
        mag_t xb;
        ulong n, k;

        arf_init(x);
        arf_init(y);
        fmpz_init(f);
        mag_init(xb);

        mag_randtest_special(xb, state, 80);
        n = n_randtest(state) % 10000;
        k = n_randtest(state) % 10000;

        mag_bin_uiui(xb, n, k);
        fmpz_bin_uiui(f, n, k);
        arf_set_fmpz(x, f);
        arf_set_mag(y, xb);

        MAG_CHECK_BITS(xb)

        if (!(arf_cmpabs(y, x) >= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        fmpz_clear(f);
        mag_clear(xb);
    }

    TEST_FUNCTION_END(state);
}
