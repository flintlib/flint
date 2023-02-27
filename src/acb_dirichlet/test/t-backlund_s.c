/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of Arb.
    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("backlund_s....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest(a, state, 1 + n_randint(state, 1000), 4);

        acb_dirichlet_backlund_s(b, a, prec1);

        if (n_randint(state, 2))
        {
            acb_dirichlet_backlund_s(c, a, prec2);
        }
        else  /* test aliasing */
        {
            arb_set(c, a);
            acb_dirichlet_backlund_s(c, c, prec2);
        }

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
