/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly_mat.h"

TEST_FUNCTION_START(fmpz_poly_mat_det, state)
{
    slong i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        fmpz_poly_t a, b, ab, c;
        slong n, bits, deg;
        float density;

        n = n_randint(state, 10);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(B, n, n);
        fmpz_poly_mat_init(C, n, n);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(ab);
        fmpz_poly_init(c);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);
        fmpz_poly_mat_randtest_sparse(B, state, deg, bits, density);
        fmpz_poly_mat_mul(C, A, B);

        fmpz_poly_mat_det(a, A);
        fmpz_poly_mat_det(b, B);
        fmpz_poly_mat_det(c, C);
        fmpz_poly_mul(ab, a, b);

        if (!fmpz_poly_equal(c, ab))
        {
            flint_printf("FAIL:\n");
            flint_printf("determinants don't agree!\n");
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            flint_printf("C:\n");
            fmpz_poly_mat_print(C, "x");
            flint_printf("det(A):\n");
            fmpz_poly_print_pretty(a, "x");
            flint_printf("\ndet(B):\n");
            fmpz_poly_print_pretty(b, "x");
            flint_printf("\ndet(C):\n");
            fmpz_poly_print_pretty(c, "x");
            flint_printf("\ndet(A)*det(B):\n");
            fmpz_poly_print_pretty(ab, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(ab);
        fmpz_poly_clear(c);

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
