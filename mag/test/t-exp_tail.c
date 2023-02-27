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

    flint_printf("exp_tail....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, t, y, z;
        fmpz_t f;
        mag_t xb, yb;
        ulong N, k;

        fmpr_init(x);
        fmpr_init(t);
        fmpr_init(y);
        fmpr_init(z);
        fmpz_init(f);
        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 6);
        mag_randtest_special(yb, state, 6);
        N = n_randint(state, 100);

        mag_exp_tail(yb, xb, N);

        mag_get_fmpr(x, xb);
        mag_get_fmpr(y, yb);

        fmpr_pow_sloppy_ui(t, x, N, MAG_BITS, FMPR_RND_DOWN);
        fmpz_fac_ui(f, N);
        fmpr_div_fmpz(t, t, f, MAG_BITS, FMPR_RND_DOWN);
        fmpr_set(z, t);

        for (k = 1; k < 50; k++)
        {
            fmpr_mul(t, t, x, MAG_BITS, FMPR_RND_DOWN);
            fmpr_div_ui(t, t, N + k, MAG_BITS, FMPR_RND_DOWN);
            fmpr_add(z, z, t, MAG_BITS, FMPR_RND_DOWN);
        }

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(fmpr_cmpabs(z, y) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("N = %wu\n\n", N);
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
            flint_abort();
        }

        mag_exp_tail(xb, xb, N);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(t);
        fmpr_clear(y);
        fmpr_clear(z);
        fmpz_clear(f);
        mag_clear(xb);
        mag_clear(yb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

