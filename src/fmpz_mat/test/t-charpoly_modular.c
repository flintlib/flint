/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_mat_charpoly_modular, state)
{
    slong rep, i;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_mat_t A, B, C, D;
        fmpz_poly_t f, g;
        slong n;

        n = n_randint(state, 8);

        flint_set_num_threads(1 + n_randint(state, 3));

        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, n, n);
        fmpz_mat_init(C, n, n);
        fmpz_mat_init(D, n, n);
        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_mat_randtest(A, state, 10);
        fmpz_mat_randtest(B, state, 10);

        fmpz_mat_mul(C, A, B);
        fmpz_mat_mul(D, B, A);

        fmpz_mat_charpoly_modular(f, C);
        fmpz_mat_charpoly_modular(g, D);

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
        slong n, bits;

        if (n_randint(state, 10) == 0)
        {
            n = n_randint(state, 4);
            bits = n_randint(state, 5000);
        }
        else
        {
            n = n_randint(state, 8);
            bits = n_randint(state, 10);
        }

        flint_set_num_threads(1 + n_randint(state, 3));

        fmpz_init(c);
        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, n, n);
        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_mat_randtest(A, state, bits);
        fmpz_mat_set(B, A);

        for (i = 0; i < 10; i++)
        {
           fmpz_randtest(c, state, 5);
           fmpz_mat_similarity(B, n_randint(state, n), c);
        }

        fmpz_mat_charpoly_modular(f, A);
        fmpz_mat_charpoly_modular(g, B);

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
