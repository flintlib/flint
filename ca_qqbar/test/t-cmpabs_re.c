/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("cmpabs_re....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        ca_qqbar_t x, y, xr, yr;
        int c1, c2;

        ca_qqbar_init(x);
        ca_qqbar_init(y);
        ca_qqbar_init(xr);
        ca_qqbar_init(yr);

        ca_qqbar_randtest(x, state, 3, 100);
        ca_qqbar_randtest(y, state, 3, 100);

        ca_qqbar_re(xr, x);
        ca_qqbar_re(yr, y);
        ca_qqbar_abs(xr, xr);
        ca_qqbar_abs(yr, yr);

        c1 = ca_qqbar_cmpabs_re(x, y);
        c2 = ca_qqbar_cmp_re(xr, yr);

        if (c1 != c2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); ca_qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); ca_qqbar_print(y); flint_printf("\n\n");
            flint_printf("%d\n\n", c1);
            flint_printf("%d\n\n", c2);
            flint_abort();
        }

        ca_qqbar_clear(x);
        ca_qqbar_clear(y);
        ca_qqbar_clear(xr);
        ca_qqbar_clear(yr);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

