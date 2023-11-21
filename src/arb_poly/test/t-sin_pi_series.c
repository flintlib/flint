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

TEST_FUNCTION_START(arb_poly_sin_pi_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n1, n2, bits1, bits2, bits3;
        arb_poly_t S, A, B, C;
        arb_t pi;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 30);
        n1 = 1 + n_randint(state, 30);
        n2 = 1 + n_randint(state, 30);

        arb_poly_init(S);
        arb_poly_init(A);
        arb_poly_init(B);
        arb_poly_init(C);
        arb_init(pi);

        arb_poly_randtest(S, state, m, bits1, 3);
        arb_poly_randtest(A, state, m, bits1, 3);
        arb_poly_randtest(B, state, m, bits1, 3);

        arb_poly_sin_pi_series(A, S, n1, bits2);

        arb_const_pi(pi, bits3);
        arb_poly_set_arb(B, pi);
        arb_poly_mul(B, S, B, bits3);
        arb_poly_sin_series(B, B, n2, bits3);

        arb_poly_set(C, A);
        arb_poly_truncate(C, FLINT_MIN(n1, n2));
        arb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!arb_poly_overlaps(B, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("S = "); arb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("A = "); arb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("B = "); arb_poly_printd(B, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_sin_pi_series(S, S, n1, bits2);

        if (!arb_poly_overlaps(A, S))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arb_poly_clear(S);
        arb_poly_clear(A);
        arb_poly_clear(B);
        arb_poly_clear(C);
        arb_clear(pi);
    }

    TEST_FUNCTION_END(state);
}
