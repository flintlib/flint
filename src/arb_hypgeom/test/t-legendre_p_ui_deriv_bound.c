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

TEST_FUNCTION_START(arb_hypgeom_legendre_p_ui_deriv_bound, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, x2sub1, y, z;
        mag_t m1, m2, exact;
        ulong n;
        slong prec;

        arb_init(x);
        arb_init(x2sub1);
        arb_init(y);
        arb_init(z);
        mag_init(m1);
        mag_init(m2);
        mag_init(exact);

        n = n_randtest(state) % 10000;
        prec = 64;

        arb_randtest(x, state, 2 * prec, 0);
        mag_zero(arb_radref(x));
        while (arf_cmpabs_2exp_si(arb_midref(x), 0) >= 0)
            arb_mul_2exp_si(x, x, -1);
        arb_mul_2exp_si(x, x, -n_randint(state, 10));

        arb_mul(x2sub1, x, x, 2 * prec);
        arb_sub_ui(x2sub1, x2sub1, 1, 2 * prec);
        arb_neg(x2sub1, x2sub1);

        arb_hypgeom_legendre_p_ui_deriv_bound(m1, m2, n, x, x2sub1);
        arb_hypgeom_legendre_p_ui(y, z, n, x, prec);
        arb_get_mag_lower(exact, z);

        if (mag_cmp(m1, exact) < 0)
        {
            flint_printf("FAIL: first derivative\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 50, 0); flint_printf("\n\n");
            flint_printf("m1 = "); mag_printd(m1, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_mul(z, z, x, prec);
        arb_mul_2exp_si(z, z, 1);
        arb_mul_ui(y, y, n, prec);
        arb_mul_ui(y, y, n + 1, prec);
        arb_sub(z, z, y, prec);
        arb_div(z, z, x2sub1, prec);
        arb_get_mag_lower(exact, z);

        if (mag_cmp(m2, exact) < 0)
        {
            flint_printf("FAIL: second derivative\n\n");
            flint_printf("n = %wu\n\n", n);
            flint_printf("x = "); arb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 50, 0); flint_printf("\n\n");
            flint_printf("m2 = "); mag_printd(m2, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(x2sub1);
        arb_clear(y);
        arb_clear(z);
        mag_clear(m1);
        mag_clear(m2);
        mag_clear(exact);
    }

    TEST_FUNCTION_END(state);
}
