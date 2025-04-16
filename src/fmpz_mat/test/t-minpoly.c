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

/* NOTE: Clang >= 18 cannot compile this on x86 systems for some reason.
 * Therefore, push this to a separate function to avoid segfaults during
 * compilation. */
FLINT_STATIC_NOINLINE void set_matrix(fmpz_mat_t am, slong n)
{
    slong ix, jx;

    for (ix = 0; ix < n/2; ix++)
        for (jx = 0; jx < n/2; jx++)
        {
            fmpz_zero(fmpz_mat_entry(am, ix + n/2, jx));
            fmpz_zero(fmpz_mat_entry(am, ix, jx + n/2));
            fmpz_set(fmpz_mat_entry(am, ix + n/2, jx + n/2), fmpz_mat_entry(am, ix, jx));
        }
}

TEST_FUNCTION_START(fmpz_mat_minpoly, state)
{
    slong n, rep;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_t c;
        fmpz_mat_t A;
        fmpz_poly_t f, g, q, r;
        slong ix;

        n = n_randint(state, 10);

        fmpz_init(c);
        fmpz_mat_init(A, n, n);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);

        if (n_randint(state, 50) == 0)
            fmpz_mat_randtest(A, state, 1000);
        else
            fmpz_mat_randtest(A, state, 10);

        set_matrix(A, n);

        for (ix = 0; ix < 10; ix++)
        {
           fmpz_randtest(c, state, 5);
           fmpz_mat_similarity(A, n_randint(state, n), c);
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
        slong ix;

        n = n_randint(state, 10);

        fmpz_init(c);
        fmpz_mat_init(A, n, n);
        fmpz_mat_init(B, n, n);
        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_mat_randtest(A, state, 10);

        set_matrix(A, n);

        fmpz_mat_set(B, A);

        for (ix = 0; ix < 10; ix++)
        {
           fmpz_randtest(c, state, 5);
           fmpz_mat_similarity(B, n_randint(state, n), c);
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

    /* special case: zero matrix used to give incorrect result: constant 1 polynomial */
    {
        fmpz_mat_t A;
        fmpz_poly_t f;

        fmpz_mat_init(A, 2, 2);
        fmpz_poly_init(f);

        fmpz_mat_zero(A);
        fmpz_mat_minpoly(f, A);

        if (fmpz_poly_length(f) != 2 || f->coeffs[0] != 0 || f->coeffs[1] != 1)
            TEST_FUNCTION_FAIL("minpoly(A) != X for zero matrix A\n"
                               "minpoly(A) = %{fmpz_poly}\n",
                               f);

        fmpz_mat_clear(A);
        fmpz_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
