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

TEST_FUNCTION_START(mag_sub, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, z2, w;
        mag_t xb, yb, zb;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(z2);
        arf_init(w);

        mag_init(xb);
        mag_init(yb);
        mag_init(zb);

        mag_randtest_special(xb, state, 100);
        mag_randtest_special(yb, state, 100);

        arf_set_mag(x, xb);
        arf_set_mag(y, yb);

        arf_sub(z, x, y, 100, ARF_RND_UP);
        if (arf_sgn(z) < 0)
            arf_zero(z);

        arf_mul_ui(z2, z, 1025, MAG_BITS, ARF_RND_UP);
        arf_mul_2exp_si(z2, z2, -10);

        mag_sub(zb, xb, yb);
        arf_set_mag(w, zb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)
        MAG_CHECK_BITS(zb)

        if (!(arf_cmpabs(z, w) <= 0 && arf_cmpabs(w, z2) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_printf("w = "); arf_print(w); flint_printf("\n\n");
            flint_abort();
        }

        if (n_randint(state, 2))
        {
            mag_sub(xb, xb, yb);

            if (!mag_equal(xb, zb))
            {
                flint_printf("FAIL (aliasing 1)\n\n");
                flint_abort();
            }
        }
        else
        {
            mag_sub(yb, xb, yb);

            if (!mag_equal(yb, zb))
            {
                flint_printf("FAIL (aliasing 2)\n\n");
                flint_abort();
            }
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(z2);
        arf_clear(w);

        mag_clear(xb);
        mag_clear(yb);
        mag_clear(zb);
    }

    TEST_FUNCTION_END(state);
}
