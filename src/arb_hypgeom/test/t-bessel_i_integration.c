/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_bessel_i_integration, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, z, r1, r2;
        slong prec1, prec2;
        int scaled;

        prec1 = 2 + n_randint(state, 80);
        prec2 = 2 + n_randint(state, 150);

        arb_init(a);
        arb_init(z);
        arb_init(r1);
        arb_init(r2);

        scaled = n_randint(state, 2);
        arb_randtest_precise(a, state, 1 + n_randint(state, 200), 1 + n_randint(state, 8));
        arb_randtest_precise(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 8));
        arb_randtest(r1, state, 1 + n_randint(state, 200), 1 + n_randint(state, 5));
        arb_randtest(r2, state, 1 + n_randint(state, 200), 1 + n_randint(state, 5));

        arb_add_ui(a, a, n_randint(state, 100), prec1 + 100);
        arb_add_ui(z, z, n_randint(state, 100), prec1 + 100);

        arb_hypgeom_bessel_i_integration(r1, a, z, scaled, prec1);

        if (arb_is_finite(r1))
        {
            if (n_randint(state, 2))
                arb_hypgeom_bessel_i_integration(r2, a, z, scaled, prec2);
            else
            {
                if (scaled)
                    arb_hypgeom_bessel_i_scaled(r2, a, z, prec2);
                else
                    arb_hypgeom_bessel_i(r2, a, z, prec2);
            }

            if (!arb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("a = "); arb_printd(a, 30); flint_printf("\n\n");
                flint_printf("z = "); arb_printd(z, 30); flint_printf("\n\n");
                flint_printf("r1 = "); arb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); arb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(a);
        arb_clear(z);
        arb_clear(r1);
        arb_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
