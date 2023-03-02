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

    flint_printf("erf_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 200 * arb_test_multiplier(); iter++)
    {
        slong m, n1, n2, bits1, bits2, bits3;
        acb_poly_t S, A, B, C;

        bits1 = 2 + n_randint(state, 200);
        bits2 = 2 + n_randint(state, 200);
        bits3 = 2 + n_randint(state, 200);

        m = 1 + n_randint(state, 10);
        n1 = 1 + n_randint(state, 10);
        n2 = 1 + n_randint(state, 10);

        acb_poly_init(S);
        acb_poly_init(A);
        acb_poly_init(B);
        acb_poly_init(C);

        acb_poly_randtest(S, state, m, bits1, 3);
        acb_poly_randtest(A, state, m, bits1, 3);
        acb_poly_randtest(B, state, m, bits1, 3);

        acb_hypgeom_erf_series(A, S, n1, bits2);
        acb_hypgeom_erf_series(B, S, n2, bits3);

        acb_poly_set(C, A);
        acb_poly_truncate(C, FLINT_MIN(n1, n2));
        acb_poly_truncate(B, FLINT_MIN(n1, n2));

        if (!acb_poly_overlaps(B, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("S = "); acb_poly_printd(S, 15); flint_printf("\n\n");
            flint_printf("A = "); acb_poly_printd(A, 15); flint_printf("\n\n");
            flint_printf("B = "); acb_poly_printd(B, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_hypgeom_erf_series(S, S, n1, bits2);

        if (!acb_poly_overlaps(A, S))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        acb_poly_clear(S);
        acb_poly_clear(A);
        acb_poly_clear(B);
        acb_poly_clear(C);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
