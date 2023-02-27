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

    flint_printf("get_mag....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b;
        mag_t m;

        arb_init(a);
        arb_init(b);
        mag_init(m);

        arb_randtest_special(a, state, 200, 1 + n_randint(state, 100));
        arb_get_mag(m, a);
        MAG_CHECK_BITS(m)

        if (arf_is_nan(arb_midref(a)))
            arf_nan(arb_midref(b));
        else
            arf_zero(arb_midref(b));
        mag_set(arb_radref(b), m);

        if (!arb_contains(b, a))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("m = "); mag_print(m); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        mag_clear(m);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

