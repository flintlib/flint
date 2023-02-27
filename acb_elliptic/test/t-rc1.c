/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_elliptic.h"
#include "acb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("rc1....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000 * arb_test_multiplier(); iter++)
    {
        acb_t x, r1, r2, a, b, c;

        acb_init(x);
        acb_init(r1);
        acb_init(r2);
        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest(x, state, 1 + n_randint(state, 300), 10);
        acb_randtest(a, state, 1 + n_randint(state, 300), 10);
        acb_add(r2, x, a, 200);
        acb_sub(r2, r2, a, 200);
        acb_neg(r2, r2);
        if (n_randint(state, 2))
            acb_swap(x, r2);

        acb_elliptic_rc1(r1, x, 2 + n_randint(state, 200));

        acb_one(a);
        acb_set_d(b, 0.5);
        acb_set_d(c, 1.5);
        acb_hypgeom_2f1(r2, a, b, c, r2, 0, 20 + n_randint(state, 200));

        if (!acb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("x = "); acb_printd(x, 30); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(r1);
        acb_clear(r2);
        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

