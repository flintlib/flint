/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_legendre_p_ui_rec, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
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
        /* close to +-1 */
        if (n_randint(state, 8) == 0)
        {
            arf_one(arb_midref(x));
            arf_mul_2exp_si(arb_midref(x), arb_midref(x), -(slong) n_randint(state, 40));
            arf_sub_ui(arb_midref(x), arb_midref(x), 1, ARF_PREC_EXACT, ARF_RND_DOWN);
            if (n_randint(state, 2))
                arf_neg(arb_midref(x), arb_midref(x));
            mag_zero(arb_radref(x));
        }
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

        /* for an exact x, the absolute error of P_n is at most about
           2^-prec1 plus the final rounding,
           and that of P_n' about n 2^-prec1 / (1 - x^2) */
        if (arb_is_exact(x) && arf_cmpabs_2exp_si(arb_midref(x), 0) < 0)
        {
            mag_t b, t;
            mag_init(b);
            mag_init(t);

            mag_set_ui_2exp_si(b, 1, 2 - prec1);
            if (mag_cmp(arb_radref(r1), b) > 0)
            {
                flint_printf("FAIL: accuracy of P_n\n\n");
                flint_printf("n = %wu, prec1 = %wd\n\n", n, prec1);
                flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
                flint_printf("r1 = "); arb_printn(r1, 50, 0); flint_printf("\n\n");
                flint_abort();
            }

            /* b = 4 n 2^-prec1 / (1 - |x|) */
            arf_get_mag(t, arb_midref(x));
            mag_one(b);
            mag_sub_lower(t, b, t);
            mag_set_ui_2exp_si(b, FLINT_MAX(n, 1), 2 - prec1);
            mag_div(b, b, t);
            if (mag_cmp(arb_radref(r3), b) > 0)
            {
                flint_printf("FAIL: accuracy of P_n'\n\n");
                flint_printf("n = %wu, prec1 = %wd\n\n", n, prec1);
                flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
                flint_printf("r3 = "); arb_printn(r3, 50, 0); flint_printf("\n\n");
                flint_abort();
            }

            mag_clear(b);
            mag_clear(t);
        }

        arb_clear(x);
        arb_clear(r1);
        arb_clear(r2);
        arb_clear(r3);
        arb_clear(r4);
    }

    TEST_FUNCTION_END(state);
}
