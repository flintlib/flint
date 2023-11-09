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
#include "fmpz_mat.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_mat_minpoly, state)
{
    slong m, n, rep, i, j;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_t c;
        fmpz_mat_t A;
        fmpz_poly_t f, g, q, r;

        m = n_randint(state, 4);
        n = m;

        fmpz_init(c);
        fmpz_mat_init(A, m, n);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);

        fmpz_mat_randtest(A, state, 10);

        for (i = 0; i < n/2; i++)
        {
           for (j = 0; j < n/2; j++)
           {
              fmpz_zero(fmpz_mat_entry(A, i + n/2, j));
              fmpz_zero(fmpz_mat_entry(A, i, j + n/2));
              fmpz_set(fmpz_mat_entry(A, i + n/2, j + n/2), fmpz_mat_entry(A, i, j));
           }
        }

        for (i = 0; i < 10; i++)
        {
           fmpz_randtest(c, state, 5);
           fmpz_mat_similarity(A, n_randint(state, m), c);
        }

        fmpz_mat_minpoly(f, A);
        fmpz_mat_charpoly(g, A);

        fmpz_poly_divrem(q, r, g, f);

        if (!fmpz_poly_is_zero(r))
        {
            flint_printf("FAIL: minpoly(A) doesn't divide charpoly(A).\n");
            flint_printf("Matrix A:\n"), fmpz_mat_print(A), flint_printf("\n");
            flint_printf("mp(A) = "), fmpz_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("cp(A) = "), fmpz_poly_print_pretty(g, "X"), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(c);
        fmpz_mat_clear(A);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
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

        for (i = 0; i < n/2; i++)
        {
           for (j = 0; j < n/2; j++)
           {
              fmpz_zero(fmpz_mat_entry(A, i + n/2, j));
              fmpz_zero(fmpz_mat_entry(A, i, j + n/2));
              fmpz_set(fmpz_mat_entry(A, i + n/2, j + n/2), fmpz_mat_entry(A, i, j));
           }
        }

        fmpz_mat_set(B, A);

        for (i = 0; i < 10; i++)
        {
           fmpz_randtest(c, state, 5);
           fmpz_mat_similarity(B, n_randint(state, m), c);
        }

        fmpz_mat_minpoly(f, A);
        fmpz_mat_minpoly(g, B);

        if (!fmpz_poly_equal(f, g))
        {
            flint_printf("FAIL: minpoly(A) != minpoly(P^{-1}AP).\n");
            flint_printf("Matrix A:\n"), fmpz_mat_print(A), flint_printf("\n");
            flint_printf("Matrix P^{-1}AP:\n"), fmpz_mat_print(B), flint_printf("\n");
            flint_printf("mp(A) = "), fmpz_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("mp(P^{-1}AP) = "), fmpz_poly_print_pretty(g, "X"), flint_printf("\n");
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
