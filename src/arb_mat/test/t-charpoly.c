/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int
main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("charpoly....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arb_mat_t A, B, C, D;
        arb_poly_t f, g;
        slong m, n;

        m = n_randint(state, 8);
        n = m;

        arb_mat_init(A, m, n);
        arb_mat_init(B, m, n);
        arb_mat_init(C, m, m);
        arb_mat_init(D, n, n);
        arb_poly_init(f);
        arb_poly_init(g);

        arb_mat_randtest(A, state, 1 + n_randint(state, 1000), 10);
        arb_mat_randtest(B, state, 1 + n_randint(state, 1000), 10);

        arb_mat_mul(C, A, B, 2 + n_randint(state, 1000));
        arb_mat_mul(D, B, A, 2 + n_randint(state, 1000));

        arb_mat_charpoly(f, C, 2 + n_randint(state, 1000));
        arb_mat_charpoly(g, D, 2 + n_randint(state, 1000));

        if (!arb_poly_overlaps(f, g))
        {
            flint_printf("FAIL: charpoly(AB) != charpoly(BA).\n");
            flint_printf("Matrix A:\n"), arb_mat_printd(A, 15), flint_printf("\n");
            flint_printf("Matrix B:\n"), arb_mat_printd(B, 15), flint_printf("\n");
            flint_printf("cp(AB) = "), arb_poly_printd(f, 15), flint_printf("\n");
            flint_printf("cp(BA) = "), arb_poly_printd(g, 15), flint_printf("\n");
            flint_abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(C);
        arb_mat_clear(D);
        arb_poly_clear(f);
        arb_poly_clear(g);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

