/*
    Copyright (C) 2014 Fredrik Johansson

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

    flint_printf("get_interval_arf....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t x, y;
        arf_t a, b;

        arb_init(x);
        arf_init(a);
        arf_init(b);
        arb_init(y);

        arb_randtest_special(x, state, 200, 100);
        arb_get_interval_arf(a, b, x, 2 + n_randint(state, 200));
        arb_set_interval_arf(y, a, b, 2 + n_randint(state, 200));

        if (!arb_contains(y, x))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("a = "); arf_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_printf("y = "); arb_print(y); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arf_clear(a);
        arf_clear(b);
        arb_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

