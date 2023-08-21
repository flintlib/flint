/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("char_get_a....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: various dots are the same */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong g = n_randint(state, 10);
        slong* n;
        ulong a, t;

        n = flint_malloc(g * sizeof(slong));
        a = n_randint(state, 1 << g);

        acb_theta_char_get_slong(n, a, g);
        t = acb_theta_char_get_a(n, g);

        if (a != t)
        {
            flint_printf("FAIL\n");
            flint_printf("a, t: %wd, %wd\n", a, t);
            flint_abort();
        }

        flint_free(n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
