/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("companion....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        acb_mat_t A;
        acb_poly_t f, g;
        slong n, prec;

        acb_poly_init(f);
        acb_poly_init(g);

        do {
            acb_poly_randtest(f, state, 1 + n_randint(state, 8), 1 + n_randint(state, 1000), 10);
        } while (acb_poly_degree(f) < 0);

        n = acb_poly_degree(f);
        prec = 2 + n_randint(state, 200);

        acb_mat_init(A, n, n);
        acb_mat_randtest(A, state, 1 + n_randint(state, 1000), 10);
        acb_mat_companion(A, f, prec);
        acb_mat_charpoly(g, A, prec);
        acb_poly_scalar_mul(g, g, acb_poly_get_coeff_ptr(f, n), prec);

        if (!acb_poly_contains(g, f))
        {
            flint_printf("FAIL\n");
            flint_printf("A:\n"), acb_mat_printd(A, 15), flint_printf("\n");
            flint_printf("f:\n"), acb_poly_printd(f, 15), flint_printf("\n");
            flint_printf("g:\n"), acb_poly_printd(g, 15), flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(A);
        acb_poly_clear(f);
        acb_poly_clear(g);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
