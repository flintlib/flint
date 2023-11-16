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

TEST_FUNCTION_START(arb_hypgeom_legendre_p_ui, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, r1, r2, r3, r4, s, t;
        ulong n;
        slong prec1, prec2, prec3;

        arb_init(x);
        arb_init(r1);
        arb_init(r2);
        arb_init(r3);
        arb_init(r4);
        arb_init(s);
        arb_init(t);

        n = n_randtest(state) % 1000;
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 500);
        prec3 = 2 + n_randint(state, 500);

        arb_randtest(x, state, 2 + n_randint(state, 1000), 3);
        arb_mul_2exp_si(x, x, -n_randint(state, 15));
        if (n_randint(state, 2))
            mag_zero(arb_radref(x));
        if (n_randint(state, 2))
            arb_add_ui(x, x, 1, prec1);

        arb_hypgeom_legendre_p_ui(r1, r2, n, x, prec1);
        if (n_randint(state, 2))
        {
            arb_hypgeom_legendre_p_ui(r3, NULL, n, x, prec2);
            arb_hypgeom_legendre_p_ui(NULL, r4, n, x, prec2);
        }
        else
        {
                arb_hypgeom_legendre_p_ui(r3, r4, n, x, prec2);
        }

        if (!arb_overlaps(r1, r3) || !arb_overlaps(r2, r4))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 50, 0); flint_printf("\n\n");
            flint_printf("r3 = "); arb_printn(r3, 50, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 50, 0); flint_printf("\n\n");
            flint_printf("r4 = "); arb_printn(r4, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (n != 0)
        {
            arb_hypgeom_legendre_p_ui(r3, NULL, n - 1, x, prec3);

            /* check (x^2-1)/n P'_n(x) = x P_n(x) - P_{n-1}(x) */

            arb_mul(s, x, x, prec1);
            arb_sub_ui(s, s, 1, prec1);
            arb_mul(s, s, r2, prec1);
            arb_div_ui(s, s, n, prec1);

            arb_mul(t, x, r1, prec1);
            arb_sub(t, t, r3, prec1);

            if (!arb_overlaps(s, t))
            {
                flint_printf("FAIL: overlap (2)\n\n");
                flint_printf("n = %wu\n\n", n);
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
        arb_clear(r4);
        arb_clear(s);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
