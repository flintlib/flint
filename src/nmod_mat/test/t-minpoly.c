/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_mat_minpoly, state)
{
    slong m, n, rep, i, j;
    ulong mod;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A;
        nmod_poly_t f, g, q, r;

        m = n_randint(state, 10);
        n = m;

        mod = n_randprime(state, 6, 0);

        nmod_mat_init(A, m, n, mod);
        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);
        nmod_poly_init(q, mod);
        nmod_poly_init(r, mod);

        nmod_mat_randtest(A, state);

        nmod_mat_minpoly(f, A);
        nmod_mat_charpoly(g, A);

        nmod_poly_divrem(q, r, g, f);

        if (!nmod_poly_is_zero(r))
        {
            flint_printf("FAIL: minpoly(A) does not divide charpoly(BA).\n");
            flint_printf("Matrix A:\n"), nmod_mat_print_pretty(A), flint_printf("\n");
            flint_printf("mp(A) = "), nmod_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("cp(A) = "), nmod_poly_print_pretty(g, "X"), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A;
        nmod_poly_t f, g;

        m = n_randint(state, 10);
        n = m;

        mod = n_randprime(state, 6, 0);

        nmod_mat_init(A, m, n, mod);
        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);

        nmod_mat_randtest(A, state);

        for (i = 0; i < n/2; i++)
        {
           for (j = 0; j < n/2; j++)
           {
              A->rows[i + n/2][j] = 0;
              A->rows[i][j + n/2] = 0;
              A->rows[i + n/2][j + n/2] = A->rows[i][j];
           }
        }

        nmod_mat_minpoly(f, A);

        for (i = 0; i < 10; i++)
           nmod_mat_similarity(A, n_randint(state, m), n_randint(state, mod));

        nmod_mat_minpoly(g, A);

        if (!nmod_poly_equal(f, g))
        {
            flint_printf("FAIL: minpoly(P^{-1}AP) != minpoly(A).\n");
            flint_printf("Matrix A:\n"), nmod_mat_print_pretty(A), flint_printf("\n");
            flint_printf("mp(A) = "), nmod_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("mp(P^{-1}AP) = "), nmod_poly_print_pretty(g, "X"), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
