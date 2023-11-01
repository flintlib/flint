/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Sebastian Pancratz

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

TEST_FUNCTION_START(nmod_mat_charpoly_danilevsky, state)
{
    slong m, n, rep, i;
    ulong mod;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C, D;
        nmod_poly_t f, g;

        m = n_randint(state, 10);
        n = m;

        mod = n_randprime(state, 6, 0);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(C, m, m, mod);
        nmod_mat_init(D, n, n, mod);
        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_mul(C, A, B);
        nmod_mat_mul(D, B, A);

        nmod_mat_charpoly_danilevsky(f, C);
        nmod_mat_charpoly_danilevsky(g, D);

        if (!nmod_poly_equal(f, g))
        {
            flint_printf("FAIL: charpoly(AB) != charpoly(BA).\n");
            flint_printf("Matrix A:\n"), nmod_mat_print_pretty(A), flint_printf("\n");
            flint_printf("Matrix B:\n"), nmod_mat_print_pretty(B), flint_printf("\n");
            flint_printf("cp(AB) = "), nmod_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("cp(BA) = "), nmod_poly_print_pretty(g, "X"), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
        nmod_poly_clear(f);
        nmod_poly_clear(g);
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

        nmod_mat_charpoly(f, A);

        for (i = 0; i < 10; i++)
           nmod_mat_similarity(A, n_randint(state, m), n_randint(state, mod));

        nmod_mat_charpoly_danilevsky(g, A);

        if (!nmod_poly_equal(f, g))
        {
            flint_printf("FAIL: charpoly(P^{-1}AP) != charpoly(A).\n");
            flint_printf("Matrix A:\n"), nmod_mat_print_pretty(A), flint_printf("\n");
            flint_printf("cp(A) = "), nmod_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("cp(P^{-1}AP) = "), nmod_poly_print_pretty(g, "X"), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
