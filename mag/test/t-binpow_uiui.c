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

    flint_printf("binpow_uiui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, y;
        mag_t b;
        ulong m, n;

        fmpr_init(x);
        fmpr_init(y);
        mag_init(b);

        m = n_randtest(state);
        n = n_randtest(state);

        fmpr_one(x);
        fmpr_div_ui(x, x, m, 128, FMPR_RND_UP);
        fmpr_add_ui(x, x, 1, 128, FMPR_RND_UP);
        fmpr_pow_sloppy_ui(x, x, n, 128, FMPR_RND_UP);

        mag_binpow_uiui(b, m, n);
        mag_get_fmpr(y, b);

        MAG_CHECK_BITS(b)

        if (!(fmpr_cmpabs(x, y) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("m = %wu\n\n", m);
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); fmpr_printd(x, 10); flint_printf("\n\n");
            flint_printf("y = "); fmpr_printd(y, 10); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        mag_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

