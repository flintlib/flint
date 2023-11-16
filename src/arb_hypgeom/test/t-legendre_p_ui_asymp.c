/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_legendre_p_ui_asymp, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, r1, r1p, r2, r2p, r3, s;
        ulong n;
        slong prec1, prec2, K, K2, K3;

        arb_init(x);
        arb_init(r1);
        arb_init(r1p);
        arb_init(r2);
        arb_init(r2p);
        arb_init(r3);
        arb_init(s);

        n = n_randtest(state) % 500;
        K = n_randint(state, 50);
        K2 = n_randint(state, 500);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 500);

        arb_randtest(x, state, 2 + n_randint(state, 1000), 0);
        arb_mul_2exp_si(x, x, -n_randint(state, 8));
        if (n_randint(state, 2))
            mag_zero(arb_radref(x));

        arb_hypgeom_legendre_p_ui_asymp(r1, r1p, n, x, K, prec1);

        arb_hypgeom_legendre_p_ui_one(r2, r2p, n, x, K2, prec2);

        if (!arb_overlaps(r1, r2) || !arb_overlaps(r1p, r2p))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu, K = %wd\n\n", n, K);
            flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 50, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 50, 0); flint_printf("\n\n");
            flint_printf("r1p = "); arb_printn(r1p, 50, 0); flint_printf("\n\n");
            flint_printf("r2p = "); arb_printn(r2p, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        n = n_randtest(state);
        if (n == UWORD_MAX) n--;
        if (n == 0) n++;
        K = n_randint(state, 50);
        K2 = n_randint(state, 50);
        K3 = n_randint(state, 50);

        /* test (2n+1) x P_n(x) - n P_(n-1)(x) = (n+1) P_(n+1)(x) */
        arb_hypgeom_legendre_p_ui_asymp(r1, NULL, n, x, K, prec2);
        arb_mul(r1, r1, x, prec2);
        arb_set_ui(r2, n);
        arb_add_ui(r2, r2, n, prec2);
        arb_add_ui(r2, r2, 1, prec2);
        arb_mul(r1, r1, r2, prec2);

        arb_hypgeom_legendre_p_ui_asymp(r2, NULL, n - 1, x, K2, prec2);
        arb_mul_ui(r2, r2, n, prec2);

        arb_hypgeom_legendre_p_ui_asymp(r3, NULL, n + 1, x, K3, prec2);
        arb_mul_ui(r3, r3, n + 1, prec2);

        arb_sub(s, r1, r2, prec2);

        if (!arb_overlaps(s, r3))
        {
            flint_printf("FAIL: overlap (2)\n\n");
            flint_printf("n = %wu, K = %wd, K2 = %wd, K3 = %wd\n\n", n, K, K2, K3);
            flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 50, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 50, 0); flint_printf("\n\n");
            flint_printf("s  = "); arb_printn(s,  50, 0); flint_printf("\n\n");
            flint_printf("r3 = "); arb_printn(r3, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(r1);
        arb_clear(r1p);
        arb_clear(r2);
        arb_clear(r2p);
        arb_clear(r3);
        arb_clear(s);
    }

    TEST_FUNCTION_END(state);
}
