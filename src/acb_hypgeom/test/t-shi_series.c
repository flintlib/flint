/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("shi_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 200 * arb_test_multiplier(); iter++)
    {
        slong m, n1, n2, n3, bits1, bits2, bits3;
        acb_poly_t S, A, B, C, T, U;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 10);
        n1 = 1 + n_randint(state, 10);
        n2 = 1 + n_randint(state, 10);
        n3 = FLINT_MIN(n1, n2);

        acb_poly_init(S);
        acb_poly_init(A);
        acb_poly_init(B);
        acb_poly_init(C);
        acb_poly_init(T);
        acb_poly_init(U);

        acb_poly_randtest(S, state, m, bits1, 3);
        acb_poly_randtest(A, state, m, bits1, 3);
        acb_poly_randtest(B, state, m, bits1, 3);

        acb_hypgeom_shi_series(A, S, n1, bits2);
        acb_hypgeom_shi_series(B, S, n2, bits3);

        acb_poly_set(C, A);
        acb_poly_truncate(C, n3);
        acb_poly_truncate(B, n3);

        /* [Shi(h(x))]' h(x) = sinh(h(x)) h'(x) */
        acb_poly_sinh_series(U, S, n3, bits2);
        acb_poly_derivative(T, S, bits2);
        acb_poly_mullow(U, U, T, FLINT_MAX(0, n3 - 1), bits2);

        acb_poly_derivative(T, A, bits2);
        acb_poly_mullow(T, T, S, FLINT_MAX(0, n3 - 1), bits2);

        if (!acb_poly_overlaps(B, C) || !acb_poly_overlaps(T, U))
        {
            flint_printf("FAIL\n\n");
            flint_printf("S = "); acb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("A = "); acb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("B = "); acb_poly_printd(B, 15); flint_printf("\n\n");
            flint_printf("T = "); acb_poly_printd(T, 15); flint_printf("\n\n");
            flint_printf("U = "); acb_poly_printd(U, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_hypgeom_shi_series(S, S, n1, bits2);

        if (!acb_poly_overlaps(A, S))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        acb_poly_clear(S);
        acb_poly_clear(A);
        acb_poly_clear(B);
        acb_poly_clear(C);
        acb_poly_clear(T);
        acb_poly_clear(U);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
