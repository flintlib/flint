/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mag.h"
#include "arb.h"

TEST_FUNCTION_START(mag_polylog_tail, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        mag_t xb, yb;
        ulong N, k, d;
        slong s, prec;
        arb_t z, t, u, r;

        mag_init(xb);
        mag_init(yb);
        arb_init(z);
        arb_init(t);
        arb_init(u);
        arb_init(r);

        mag_randtest_special(xb, state, 6);
        mag_randtest_special(yb, state, 6);
        N = n_randint(state, 100);
        d = n_randint(state, 100);
        s = n_randint(state, 100) - 50;
        prec = 4 * MAG_BITS;

        mag_polylog_tail(yb, xb, s, d, N);

        arb_zero(z);
        arf_set_mag(arb_midref(z), xb);
        arb_zero(r);

        for (k = N; k < N + 100; k++)
        {
            arb_pow_ui(t, z, k, prec);
            arb_log_ui(u, k, prec);
            arb_pow_ui(u, u, d, prec);
            arb_mul(t, t, u, prec);

            arb_set_ui(u, k);
            if (s >= 0)
            {
                arb_pow_ui(u, u, s, prec);
                arb_div(t, t, u, prec);
            }
            else
            {
                arb_pow_ui(u, u, -s, prec);
                arb_mul(t, t, u, prec);
            }

            arb_add(r, r, t, prec);
        }

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        arb_zero(z);
        mag_set(arb_radref(z), yb);
        if (!arb_is_finite(z))
            arb_indeterminate(z);

        if (!arb_contains(z, r))
        {
            flint_printf("FAIL\n\n");
            flint_printf("N = %wu\n\n", N);
            flint_printf("d = %wu\n\n", d);
            flint_printf("s = %wd\n\n", s);
            flint_printf("xb = "); mag_printd(xb, 15); flint_printf("\n\n");
            flint_printf("yb = "); mag_printd(yb, 15); flint_printf("\n\n");
            flint_printf("z = "); arb_printd(z, 15); flint_printf("\n\n");
            flint_printf("r = "); arb_printd(r, 15); flint_printf("\n\n");
            flint_abort();
        }

        mag_polylog_tail(xb, xb, s, d, N);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        mag_clear(xb);
        mag_clear(yb);
        arb_clear(z);
        arb_clear(t);
        arb_clear(u);
        arb_clear(r);
    }

    TEST_FUNCTION_END(state);
}
