/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mag.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("root....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, y, z, z2;
        mag_t xb, yb;
        ulong n;

        fmpr_init(x);
        fmpr_init(y);
        fmpr_init(z);
        fmpr_init(z2);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 1 + n_randint(state, 200));
        mag_randtest_special(yb, state, 1 + n_randint(state, 200));
        n = n_randint(state, 30);

        mag_root(yb, xb, n);

        mag_get_fmpr(x, xb);
        mag_get_fmpr(y, yb);

        fmpr_root(z, x, n, MAG_BITS, FMPR_RND_UP);
        fmpr_mul_ui(z2, z, 1025, MAG_BITS, FMPR_RND_UP);
        fmpr_mul_2exp_si(z2, z2, -10);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(fmpr_cmpabs(z, y) <= 0 && fmpr_cmpabs(y, z2) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_printf("z2 = "); fmpr_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        mag_root(xb, xb, n);

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
