/*
    Copyright (C) 2017 Fredrik Johansson

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

    flint_printf("set_interval_mag....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        arb_t x;
        mag_t u, v;
        arf_t a, b;

        arb_init(x);
        mag_init(u);
        mag_init(v);
        arf_init(a);
        arf_init(b);

        mag_randtest_special(u, state, 1 + n_randint(state, 100));
        mag_randtest_special(v, state, 1 + n_randint(state, 100));
        if (mag_cmp(u, v) > 0)
            mag_swap(u, v);
        arf_set_mag(a, u);
        arf_set_mag(b, v);

        arb_set_interval_mag(x, u, v, 2 + n_randint(state, 200));

        if (!arb_contains_arf(x, a) || !arb_contains_arf(x, b))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("x = "); arb_print(x); flint_printf("\n\n");
            flint_printf("a = "); arf_print(a); flint_printf("\n\n");
            flint_printf("b = "); arf_print(b); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arf_clear(a);
        arf_clear(b);
        mag_clear(u);
        mag_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

