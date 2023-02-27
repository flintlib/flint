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

    flint_printf("fac_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, y;
        fmpz_t f;
        mag_t xb;
        ulong n;

        fmpr_init(x);
        fmpr_init(y);
        fmpz_init(f);
        mag_init(xb);

        mag_randtest_special(xb, state, 80);
        n = n_randtest(state) % 2000;

        mag_fac_ui(xb, n);
        fmpz_fac_ui(f, n);
        fmpr_set_fmpz(x, f);
        mag_get_fmpr(y, xb);

        MAG_CHECK_BITS(xb)

        if (!(fmpr_cmpabs(y, x) >= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);
        fmpz_clear(f);
        mag_clear(xb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

