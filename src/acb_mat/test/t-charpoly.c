/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_mat.h"

TEST_FUNCTION_START(acb_mat_charpoly, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, B, C, D;
        acb_poly_t f, g;
        slong m, n;

        m = n_randint(state, 8);
        n = m;

        acb_mat_init(A, m, n);
        acb_mat_init(B, m, n);
        acb_mat_init(C, m, m);
        acb_mat_init(D, n, n);
        acb_poly_init(f);
        acb_poly_init(g);

        acb_mat_randtest(A, state, 1 + n_randint(state, 1000), 10);
        acb_mat_randtest(B, state, 1 + n_randint(state, 1000), 10);

        acb_mat_mul(C, A, B, 2 + n_randint(state, 1000));
        acb_mat_mul(D, B, A, 2 + n_randint(state, 1000));

        acb_mat_charpoly(f, C, 2 + n_randint(state, 1000));
        acb_mat_charpoly(g, D, 2 + n_randint(state, 1000));

        if (!acb_poly_overlaps(f, g))
        {
            flint_printf("FAIL: charpoly(AB) != charpoly(BA).\n");
            flint_printf("Matrix A:\n"), acb_mat_printd(A, 15), flint_printf("\n");
            flint_printf("Matrix B:\n"), acb_mat_printd(B, 15), flint_printf("\n");
            flint_printf("cp(AB) = "), acb_poly_printd(f, 15), flint_printf("\n");
            flint_printf("cp(BA) = "), acb_poly_printd(g, 15), flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(A);
        acb_mat_clear(B);
        acb_mat_clear(C);
        acb_mat_clear(D);
        acb_poly_clear(f);
        acb_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
