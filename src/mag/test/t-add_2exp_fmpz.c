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

TEST_FUNCTION_START(mag_add_2exp_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, z2, w;
        mag_t xb, zb;
        fmpz_t e;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(z2);
        arf_init(w);

        fmpz_init(e);

        mag_init(xb);
        mag_init(zb);

        arf_randtest(x, state, MAG_BITS, 100);
        arf_abs(x, x);

        fmpz_randtest(e, state, 100);
        arf_one(y);
        arf_mul_2exp_fmpz(y, y, e);

        arf_add(z, x, y, MAG_BITS + 10, ARF_RND_DOWN);
        arf_mul_ui(z2, z, 1025, MAG_BITS, ARF_RND_UP);
        arf_mul_2exp_si(z2, z2, -10);

        arf_get_mag(xb, x);

        mag_add_2exp_fmpz(zb, xb, e);
        arf_set_mag(w, zb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(zb)

        if (!(arf_cmpabs(z, w) <= 0 && arf_cmpabs(w, z2) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); arf_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); arf_printd(z, 15); flint_printf("\n\n");
            flint_printf("w = "); arf_printd(w, 15); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(z2);
        arf_clear(w);

        mag_clear(xb);
        mag_clear(zb);

        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
