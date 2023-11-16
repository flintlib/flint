/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "arf.h"
#include "mag.h"

void
arf_sinh(arf_t y, const arf_t x, slong prec, arf_rnd_t rnd)
{
    _arf_call_mpfr_func(y, NULL, (int (*)(void)) mpfr_sinh, x, NULL, prec, rnd);
}

TEST_FUNCTION_START(mag_sinh, state)
{
    slong iter;

    /* test sinh */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, z2;
        mag_t xb, yb;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(z2);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 0);
        mag_mul_2exp_si(xb, xb, -100 + n_randint(state,120));
        mag_randtest_special(yb, state, 25);

        mag_sinh(yb, xb);

        arf_set_mag(x, xb);
        arf_set_mag(y, yb);

        arf_sinh(z, x, MAG_BITS, ARF_RND_UP);
        arf_mul_ui(z2, z, 1025, MAG_BITS, ARF_RND_UP);
        arf_mul_2exp_si(z2, z2, -10);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(arf_cmpabs(z, y) <= 0 && arf_cmpabs(y, z2) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); arf_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        mag_sinh(xb, xb);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(z2);

        mag_clear(xb);
        mag_clear(yb);
    }

    /* test sinh_lower */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, z2;
        mag_t xb, yb;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(z2);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 0);
        mag_mul_2exp_si(xb, xb, -100 + n_randint(state,120));
        mag_randtest_special(yb, state, 25);

        mag_sinh_lower(yb, xb);

        arf_set_mag(x, xb);
        arf_set_mag(y, yb);

        arf_sinh(z, x, MAG_BITS, ARF_RND_DOWN);
        arf_mul_ui(z2, z, 1023, MAG_BITS, ARF_RND_DOWN);
        arf_mul_2exp_si(z2, z2, -10);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(arf_cmpabs(z2, y) <= 0 && arf_cmpabs(y, z) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); arf_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        mag_sinh_lower(xb, xb);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(z2);

        mag_clear(xb);
        mag_clear(yb);
    }

    TEST_FUNCTION_END(state);
}
