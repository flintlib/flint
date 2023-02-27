/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("union....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t x, y, z;
        slong prec;
        int alias;

        arb_init(x);
        arb_init(y);
        arb_init(z);

        arb_randtest_special(x, state, 200, 10);
        arb_randtest_special(y, state, 200, 10);
        arb_randtest_special(z, state, 200, 10);

        prec = 2 + n_randint(state, 200);

        arb_union(z, x, y, prec);

        if (!arb_contains(z, x) || !arb_contains(z, y))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_abort();
        }

        if (n_randint(state, 2))
        {
            arb_union(x, x, y, prec);
            alias = arb_equal(x, z);
        }
        else
        {
            arb_union(y, x, y, prec);
            alias = arb_equal(y, z);
        }

        if (!alias)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

