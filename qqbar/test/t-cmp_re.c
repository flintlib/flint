/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("cmp_re....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x, y, xr, yr, t, u;
        int c1, c2;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(xr);
        qqbar_init(yr);
        qqbar_init(t);
        qqbar_init(u);

        qqbar_randtest(x, state, 3, 100);
        if (n_randint(state, 2))
            qqbar_set(y, x);
        else
            qqbar_randtest(y, state, 3, 10);
        if (n_randint(state, 2))
            qqbar_conj(y, y);

        qqbar_randtest_real(t, state, 1, 100);
        qqbar_i(u);
        qqbar_mul(t, t, u);
        qqbar_add(x, x, t);

        qqbar_re(xr, x);
        qqbar_re(yr, y);
        qqbar_sub(t, xr, yr);

        c1 = qqbar_cmp_re(x, y);
        c2 = qqbar_sgn_re(t);

        if (c1 != c2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("%d\n\n", c1);
            flint_printf("%d\n\n", c2);
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(xr);
        qqbar_clear(yr);
        qqbar_clear(t);
        qqbar_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

