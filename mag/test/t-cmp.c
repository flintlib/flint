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

    flint_printf("cmp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        fmpr_t x, y;
        mag_t xb, yb;
        int c1, c2;

        fmpr_init(x);
        fmpr_init(y);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 100);
        mag_randtest_special(yb, state, 100);

        mag_get_fmpr(x, xb);
        mag_get_fmpr(y, yb);

        c1 = fmpr_cmp(x, y);
        c2 = mag_cmp(xb, yb);

        if (c1 != c2)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpr_print(y); flint_printf("\n\n");
            flint_printf("xb = "); mag_print(xb); flint_printf("\n\n");
            flint_printf("yb = "); mag_print(yb); flint_printf("\n\n");
            flint_abort();
        }

        fmpr_clear(x);
        fmpr_clear(y);

        mag_clear(xb);
        mag_clear(yb);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

