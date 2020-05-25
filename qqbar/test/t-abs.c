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

    flint_printf("abs....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        qqbar_t x, y, z, t, u;

        qqbar_init(x);
        qqbar_init(y);
        qqbar_init(z);
        qqbar_init(t);
        qqbar_init(u);

        qqbar_randtest(x, state, 3, 10);
        qqbar_randtest(y, state, 3, 10);

        qqbar_abs(z, x);
        qqbar_abs(t, y);
        qqbar_mul(t, t, z);

        qqbar_mul(u, x, y);
        qqbar_abs(u, u);

        if (!qqbar_equal(t, u) || !qqbar_is_real(t))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_printf("z = "); qqbar_print(z); flint_printf("\n\n");
            flint_printf("t = "); qqbar_print(t); flint_printf("\n\n");
            flint_printf("u = "); qqbar_print(t); flint_printf("\n\n");
            flint_abort();
        }

        qqbar_clear(x);
        qqbar_clear(y);
        qqbar_clear(z);
        qqbar_clear(t);
        qqbar_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

