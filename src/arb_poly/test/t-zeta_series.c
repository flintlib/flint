/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_zeta_series, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n1, n2, bits1, bits2, bits3;
        int deflate;
        arb_poly_t S, A, B, C, D, E, F;
        arb_t a, a1;

        bits1 = 2 + n_randint(state, 300);
        bits2 = 2 + n_randint(state, 300);
        bits3 = 2 + n_randint(state, 300);

        m = 1 + n_randint(state, 30);
        n1 = 1 + n_randint(state, 30);
        n2 = 1 + n_randint(state, 30);

        arb_poly_init(S);
        arb_poly_init(A);
        arb_poly_init(B);
        arb_poly_init(C);
        arb_poly_init(D);
        arb_poly_init(E);
        arb_poly_init(F);
        arb_init(a);
        arb_init(a1);

        deflate = n_randint(state, 2);

        arb_poly_randtest(S, state, m, bits1, 3);
        arb_randtest_precise(a, state, bits1, 3);
        arb_poly_set_coeff_arb(S, 0, a);

        if (n_randint(state, 2))
            arb_randtest(a, state, bits1, 3);
        else
            arb_one(a);

        arb_poly_zeta_series(A, S, a, deflate, n1, bits2);
        arb_poly_zeta_series(B, S, a, deflate, n2, bits3);

        arb_poly_set(C, A);
        arb_poly_truncate(C, FLINT_MIN(n1, n2));
        arb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!arb_poly_overlaps(B, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("S = "); arb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("A = "); arb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("B = "); arb_poly_printd(B, 15); flint_printf("\n\n");
            flint_abort();
        }

        /* check zeta(s,a) = zeta(s,a+1) + a^(-s) */
        arb_poly_set_arb(D, a);
        arb_poly_log_series(D, D, n1, bits2);
        arb_poly_mullow(D, D, S, n1, bits2);
        arb_poly_neg(D, D);
        arb_poly_exp_series(D, D, n1, bits2);

        arb_add_ui(a1, a, 1, bits2);
        arb_poly_zeta_series(E, S, a1, deflate, n1, bits2);
        arb_poly_add(E, E, D, bits2);

        if (!arb_poly_overlaps(A, E))
        {
            flint_printf("FAIL (functional equation)\n\n");
            flint_printf("S = "); arb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("A = "); arb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("E = "); arb_poly_printd(A, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_zeta_series(S, S, a, deflate, n1, bits2);
        if (!arb_poly_overlaps(A, S))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arb_poly_clear(S);
        arb_poly_clear(A);
        arb_poly_clear(B);
        arb_poly_clear(C);
        arb_poly_clear(D);
        arb_poly_clear(E);
        arb_poly_clear(F);
        arb_clear(a);
        arb_clear(a1);
    }

    TEST_FUNCTION_END(state);
}
