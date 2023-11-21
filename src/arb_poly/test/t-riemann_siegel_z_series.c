/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_riemann_siegel_z_series, state)
{
    slong iter;

    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n1, n2, rbits1, rbits2, rbits3;
        arb_poly_t a, b, c, d;

        rbits1 = 2 + n_randint(state, 150);
        rbits2 = 2 + n_randint(state, 150);
        rbits3 = 2 + n_randint(state, 150);

        m = 1 + n_randint(state, 15);
        n1 = 1 + n_randint(state, 15);
        n2 = 1 + n_randint(state, 15);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_init(d);

        arb_poly_randtest(a, state, m, rbits1, 10);
        if (a->length != 0)
            arb_randtest_precise(a->coeffs, state, rbits1, 1);
        arb_poly_randtest(b, state, m, rbits1, 10);
        arb_poly_randtest(c, state, m, rbits1, 10);

        arb_poly_riemann_siegel_z_series(b, a, n1, rbits2);
        arb_poly_riemann_siegel_z_series(c, a, n2, rbits3);

        arb_poly_set(d, b);
        arb_poly_truncate(d, FLINT_MIN(n1, n2));
        arb_poly_truncate(c, FLINT_MIN(n1, n2));

        if (!arb_poly_overlaps(c, d))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n1 = %wd, n2 = %wd, bits2 = %wd, bits3 = %wd\n", n1, n2, rbits2, rbits3);

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        /* check Z(a) = Z(-a) */
        arb_poly_neg(d, a);
        arb_poly_riemann_siegel_z_series(c, d, n1, rbits2);

        if (!arb_poly_overlaps(b, c))
        {
            flint_printf("FAIL (symmetry, n1 = %wd)\n\n", n1);

            flint_printf("a = "); arb_poly_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_poly_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_poly_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); arb_poly_printd(d, 15); flint_printf("\n\n");

            flint_abort();
        }

        arb_poly_riemann_siegel_z_series(a, a, n1, rbits2);
        if (!arb_poly_overlaps(a, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
        arb_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}
