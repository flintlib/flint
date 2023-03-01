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

    flint_printf("legendre_p_ui_one....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000 * arb_test_multiplier(); iter++)
    {
        arb_t x, r1, r2, r3, s, t;
        ulong n;
        slong prec1, prec2, prec3, K, K2;

        arb_init(x);
        arb_init(r1);
        arb_init(r2);
        arb_init(r3);
        arb_init(s);
        arb_init(t);

        n = n_randtest(state) % 500;
        K = n_randint(state, 300);
        K2 = n_randint(state, 300);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 500);
        prec3 = 2 + n_randint(state, 500);

        arb_randtest(x, state, 2 + n_randint(state, 1000), 1);
        arb_mul_2exp_si(x, x, -n_randint(state, 15));
        if (n_randint(state, 2))
            mag_zero(arb_radref(x));
        if (n_randint(state, 2))
            arb_add_ui(x, x, 1, prec1);

        arb_set_ui(r1, n);
        arb_set_ui(r2, 0);
        /* todo: replace with a different implementation later */
        arb_hypgeom_legendre_p(r2, r1, r2, x, 0, prec1);

        arb_hypgeom_legendre_p_ui_one(r1, r3, n, x, K, prec2);

        if (!arb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu, K = %wd\n\n", n, K);
            flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 50, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (n != 0)
        {
            arb_hypgeom_legendre_p_ui_one(r2, NULL, n - 1, x, K2, prec3);

            /* check (x^2-1)/n P'_n(x) = x P_n(x) - P_{n-1}(x) */

            arb_mul(s, x, x, prec1);
            arb_sub_ui(s, s, 1, prec1);
            arb_mul(s, s, r3, prec1);
            arb_div_ui(s, s, n, prec1);

            arb_mul(t, x, r1, prec1);
            arb_sub(t, t, r2, prec1);

            if (!arb_overlaps(s, t))
            {
                flint_printf("FAIL: overlap (2)\n\n");
                flint_printf("n = %wu, K = %wd, K2 = %wd\n\n", n, K, K2);
                flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
                flint_printf("r1 = "); arb_printn(r1, 50, 0); flint_printf("\n\n");
                flint_printf("r2 = "); arb_printn(r2, 50, 0); flint_printf("\n\n");
                flint_printf("r3 = "); arb_printn(r3, 50, 0); flint_printf("\n\n");
                flint_printf("s  = "); arb_printn(s,  50, 0); flint_printf("\n\n");
                flint_printf("t  = "); arb_printn(t,  50, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(x);
        arb_clear(r1);
        arb_clear(r2);
        arb_clear(r3);
        arb_clear(s);
        arb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

