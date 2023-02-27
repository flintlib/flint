/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mag.h"

slong
fmpr_sinh(fmpr_t y, const fmpr_t x, slong prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        fmpr_set(y, x);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        slong r;
        CALL_MPFR_FUNC(r, mpfr_sinh, y, x, prec, rnd);
        return r;
    }
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("sinh_lower....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, y, z, z2;
        mag_t xb, yb;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(z2);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 0);
        mag_mul_2exp_si(xb, xb, -100 + n_randint(state,120));
        mag_randtest_special(yb, state, 25);

        mag_sinh_lower(yb, xb);

        mag_get_fmpr(x, xb);
        mag_get_fmpr(y, yb);

        fmpr_sinh(z, x, MAG_BITS, FMPR_RND_DOWN);
        fmpr_mul_ui(z2, z, 1023, MAG_BITS, FMPR_RND_DOWN);
        fmpr_mul_2exp_si(z2, z2, -10);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(fmpr_cmpabs(z2, y) <= 0 && fmpr_cmpabs(y, z) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); fmpr_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        mag_sinh_lower(xb, xb);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpr_clear(z2);

        mag_clear(xb);
        mag_clear(yb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

