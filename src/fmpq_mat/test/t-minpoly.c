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
#include "fmpq.h"
#include "fmpq_poly.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_minpoly, state)
{
    slong m, n, rep, i, j;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpq_mat_t A;
        fmpq_poly_t f, g, q, r;

        m = n_randint(state, 4);
        n = m;

        fmpq_mat_init(A, m, n);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(q);
        fmpq_poly_init(r);

        fmpq_mat_randtest(A, state, 10);

        fmpq_mat_charpoly(f, A);
        fmpq_mat_minpoly(g, A);

        fmpq_poly_divrem(q, r, f, g);

        if (!fmpq_poly_is_zero(r))
        {
            flint_printf("FAIL: minpoly(A) doesn't divide charpoly(A).\n");
            flint_printf("Matrix A:\n"), fmpq_mat_print(A), flint_printf("\n");
            flint_printf("cp(A) = "), fmpq_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("mp(A) = "), fmpq_poly_print_pretty(g, "X"), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(q);
        fmpq_poly_clear(r);
    }

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpq_mat_t A, B;
        fmpq_poly_t f, g;
        fmpq_t d;

        m = n_randint(state, 4);
        n = m;

        fmpq_init(d);
        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpq_poly_init(f);
        fmpq_poly_init(g);

        fmpq_mat_randtest(A, state, 10);

        for (i = 0; i < n/2; i++)
        {
           for (j = 0; j < n/2; j++)
           {
              fmpq_zero(fmpq_mat_entry(A, i, j + n/2));
              fmpq_zero(fmpq_mat_entry(A, i + n/2, j));
           }
        }

        for (i = 0; i < n/2; i++)
        {
           for (j = 0; j < n/2; j++)
           {
              fmpq_set(fmpq_mat_entry(A, i + n/2, j + n/2), fmpq_mat_entry(A, i, j));
           }
        }

        fmpq_mat_set(B, A);

        fmpq_mat_minpoly(g, A);

        for (i = 0; i < n; i++)
        {
           fmpq_set_si(d, n_randint(state, 6) - 3, 1);
           fmpq_mat_similarity(B, n_randint(state, n), d);
        }

        fmpq_mat_minpoly(f, B);

        if (!fmpq_poly_equal(f, g))
        {
            flint_printf("FAIL: minpoly(P^{-1}AP) != minpoly(A).\n");
            flint_printf("Matrix A:\n"), fmpq_mat_print(A), flint_printf("\n");
            flint_printf("Matrix P^{-1}AP:\n"), fmpq_mat_print(B), flint_printf("\n");
            flint_printf("mp(P^{-1}AP) = "), fmpq_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("mp(A) = "), fmpq_poly_print_pretty(g, "X"), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(d);
        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
