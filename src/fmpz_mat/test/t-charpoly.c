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
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_charpoly, state)
{
    slong m, n, rep, i;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_mat_t A, B, C, D;
        fmpz_poly_t f, g;

        m = n_randint(state, 4);
        n = m;

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(C, m, m);
        fmpz_mat_init(D, n, n);
        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_mat_randtest(A, state, 10);
        fmpz_mat_randtest(B, state, 10);

        fmpz_mat_mul(C, A, B);
        fmpz_mat_mul(D, B, A);

        fmpz_mat_charpoly(f, C);
        fmpz_mat_charpoly(g, D);

        if (!fmpz_poly_equal(f, g))
        {
            flint_printf("FAIL: charpoly(AB) != charpoly(BA).\n");
            flint_printf("Matrix A:\n"), fmpz_mat_print(A), flint_printf("\n");
            flint_printf("Matrix B:\n"), fmpz_mat_print(B), flint_printf("\n");
            flint_printf("cp(AB) = "), fmpz_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("cp(BA) = "), fmpz_poly_print_pretty(g, "X"), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_t c;
        fmpz_mat_t A, B;
        fmpz_poly_t f, g;

        m = n_randint(state, 4);
        n = m;

        fmpz_init(c);
        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, m, n);
        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_mat_randtest(A, state, 10);
        fmpz_mat_set(B, A);

        for (i = 0; i < 10; i++)
        {
           fmpz_randtest(c, state, 5);
           fmpz_mat_similarity(B, n_randint(state, m), c);
        }

        fmpz_mat_charpoly(f, A);
        fmpz_mat_charpoly(g, B);

        if (!fmpz_poly_equal(f, g))
        {
            flint_printf("FAIL: charpoly(A) != charpoly(P^{-1}AP).\n");
            flint_printf("Matrix A:\n"), fmpz_mat_print(A), flint_printf("\n");
            flint_printf("Matrix P^{-1}AP:\n"), fmpz_mat_print(B), flint_printf("\n");
            flint_printf("cp(A) = "), fmpz_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("cp(P^{-1}AP) = "), fmpz_poly_print_pretty(g, "X"), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(c);
        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
