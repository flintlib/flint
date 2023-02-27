/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("legendre_p_ui_rec....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000 * arb_test_multiplier(); iter++)
    {
        arb_t x, r1, r2, r3, r4;
        ulong n;
        slong prec1, prec2;

        arb_init(x);
        arb_init(r1);
        arb_init(r2);
        arb_init(r3);
        arb_init(r4);

        n = n_randtest(state) % 500;
        prec1 = 2 + n_randint(state, 500);
        prec2 = 2 + n_randint(state, 500);

        arb_randtest(x, state, 2 + n_randint(state, 1000), 1);
        arb_mul_2exp_si(x, x, -n_randint(state, 15));
        if (n_randint(state, 2))
            mag_zero(arb_radref(x));

        arb_hypgeom_legendre_p_ui_rec(r1, r3, n, x, prec1);
        arb_hypgeom_legendre_p_ui_zero(r2, r4, n, x, n + 1, prec2);

        if (!arb_overlaps(r1, r2) || !arb_overlaps(r3, r4))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 50, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 50, 0); flint_printf("\n\n");
            flint_printf("r3 = "); arb_printn(r3, 50, 0); flint_printf("\n\n");
            flint_printf("r4 = "); arb_printn(r4, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(r1);
        arb_clear(r2);
        arb_clear(r3);
        arb_clear(r4);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

