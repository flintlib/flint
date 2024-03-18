/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpfr.h"
#include "arf.h"
#include "mag.h"

void
arf_expm1(arf_t y, const arf_t x, slong prec, arf_rnd_t rnd)
{
    _arf_call_mpfr_func(y, NULL, (int (*)(void)) mpfr_expm1, x, NULL, prec, rnd);
}

TEST_FUNCTION_START(mag_expm1, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
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
        mag_randtest_special(yb, state, 0);

        mag_mul_2exp_si(xb, xb, -100 + n_randint(state,120));

        mag_expm1(yb, xb);

        arf_set_mag(x, xb);
        arf_set_mag(y, yb);

        arf_expm1(z, x, MAG_BITS, ARF_RND_UP);

        if (arf_cmpabs_ui(x, 1000) < 0)
        {
            arf_mul_ui(z2, z, 1025, MAG_BITS, ARF_RND_UP);
            arf_mul_2exp_si(z2, z2, -10);
        }
        else
        {
            arf_mul_2exp_si(z2, z, 2);
        }

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(arf_cmpabs(z, y) <= 0 && arf_cmpabs(y, z2) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); arf_printd(x, 15); flint_printf("\n\n");
            flint_printf("y = "); arf_printd(y, 15); flint_printf("\n\n");
            flint_printf("z = "); arf_printd(z, 15); flint_printf("\n\n");
            flint_printf("z2 = "); arf_printd(z2, 15); flint_printf("\n\n");
            flint_abort();
        }

        mag_expm1(xb, xb);

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
