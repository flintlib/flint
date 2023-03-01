/*
    Copyright (C) 2015 Fredrik Johansson

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

    flint_printf("contains_int....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t a;
        slong c;
        int r, ok;

        arb_init(a);
        arb_randtest_special(a, state, 1 + n_randint(state, 500), 2);

        r = arb_contains_int(a);

        ok = !r;
        for (c = 0; c < 10; c++)
        {
            if (arb_contains_si(a, c) || arb_contains_si(a, -c))
            {
                ok = !ok;
                break;
            }
        }

        if (!ok)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "); arb_printd(a, 30); flint_printf("\n\n");
            flint_printf("r = %d\n\n", r);
            flint_abort();
        }

        arb_clear(a);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
