/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_fresnel_series, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n1, n2, bits1, bits2, bits3;
        acb_poly_t X, S1, S2, C1, C2, T;
        acb_t c;
        int normalized;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 10);
        n1 = 1 + n_randint(state, 10);
        n2 = 1 + n_randint(state, 10);
        normalized = n_randint(state, 2);

        acb_poly_init(X);
        acb_poly_init(S1);
        acb_poly_init(S2);
        acb_poly_init(C1);
        acb_poly_init(C2);
        acb_poly_init(T);
        acb_init(c);

        acb_poly_randtest(X, state, m, bits1, 3);
        acb_poly_randtest(S1, state, 10, bits1, 3);
        acb_poly_randtest(S2, state, 10, bits1, 3);
        acb_poly_randtest(C1, state, 10, bits1, 3);
        acb_poly_randtest(C2, state, 10, bits1, 3);

        acb_hypgeom_fresnel_series(S1, C1, X, normalized, n1, bits2);

        switch (n_randint(state, 8))
        {
            case 0:
                acb_hypgeom_fresnel_series(S2, C2, X, normalized, n2, bits3);
                break;
            case 1:
                acb_hypgeom_fresnel_series(S2, NULL, X, normalized, n2, bits3);
                acb_hypgeom_fresnel_series(NULL, C2, X, normalized, n2, bits3);
                break;
            case 2:
                acb_poly_set(S2, X);
                acb_hypgeom_fresnel_series(NULL, C2, S2, normalized, n2, bits3);
                acb_hypgeom_fresnel_series(S2, NULL, S2, normalized, n2, bits3);
                break;
            case 3:
                acb_poly_set(C2, X);
                acb_hypgeom_fresnel_series(S2, NULL, C2, normalized, n2, bits3);
                acb_hypgeom_fresnel_series(NULL, C2, C2, normalized, n2, bits3);
                break;
            case 4:
                acb_poly_set(S2, X);
                acb_hypgeom_fresnel_series(S2, C2, S2, normalized, n2, bits3);
                break;
            case 5:
                acb_poly_set(C2, X);
                acb_hypgeom_fresnel_series(S2, C2, C2, normalized, n2, bits3);
                break;
            default:
                acb_const_pi(c, bits3);
                acb_mul_2exp_si(c, c, -1);
                acb_sqrt(c, c, bits3);

                if (normalized == 0)
                    acb_inv(c, c, bits3);

                acb_poly_set_acb(T, c);
                acb_poly_mul(T, T, X, bits3);

                acb_hypgeom_fresnel_series(S2, C2, T, !normalized, n2, bits3);

                acb_inv(c, c, bits3);
                acb_poly_set_acb(T, c);
                acb_poly_mul(S2, S2, T, bits3);
                acb_poly_mul(C2, C2, T, bits3);
        }

        acb_poly_truncate(S1, FLINT_MIN(n1, n2));
        acb_poly_truncate(C1, FLINT_MIN(n1, n2));
        acb_poly_truncate(S2, FLINT_MIN(n1, n2));
        acb_poly_truncate(C2, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(S1, S2) || !acb_poly_overlaps(C1, C2))
        {
            flint_printf("FAIL (n1 = %wd, n2 = %wd, norm = %d)\n\n", n1, n2, normalized);
            flint_printf("X = "); acb_poly_printd(X, 15); flint_printf("\n\n");
            flint_printf("S1 = "); acb_poly_printd(S1, 15); flint_printf("\n\n");
            flint_printf("S2 = "); acb_poly_printd(S2, 15); flint_printf("\n\n");
            flint_printf("C1 = "); acb_poly_printd(C1, 15); flint_printf("\n\n");
            flint_printf("C2 = "); acb_poly_printd(C2, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(X);
        acb_poly_clear(S1);
        acb_poly_clear(S2);
        acb_poly_clear(C1);
        acb_poly_clear(C2);
        acb_poly_clear(T);
        acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
