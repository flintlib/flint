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

TEST_FUNCTION_START(fmpz_poly_mat_solve_fflu, state)
{
    slong i;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, X, B, AX, Bden;
        fmpz_poly_t den, det;
        slong n, m, bits, deg;
        float density;
        int solved;

        n = n_randint(state, 15);
        m = n_randint(state, 5);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(B, n, m);
        fmpz_poly_mat_init(X, n, m);
        fmpz_poly_mat_init(AX, n, m);
        fmpz_poly_mat_init(Bden, n, m);
        fmpz_poly_init(den);
        fmpz_poly_init(det);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);
        fmpz_poly_mat_randtest_sparse(B, state, deg, bits, density);

        solved = fmpz_poly_mat_solve_fflu(X, den, A, B);
        fmpz_poly_mat_det_interpolate(det, A);

        if (m == 0 || n == 0)
        {
            if (solved == 0)
            {
                flint_printf("FAIL: expected empty system to pass\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
	    if (!fmpz_poly_equal(den, det))
            {
                fmpz_poly_neg(det, det);
                flint_printf("FAIL: den != +/- det(A)\n");
                flint_printf("den:\n"); fmpz_poly_print_pretty(den, "x");
                flint_printf("\n\n");
                flint_printf("det:\n"); fmpz_poly_print_pretty(det, "x");
                flint_printf("\n\n");
                flint_printf("A:\n");
                fmpz_poly_mat_print(A, "x");
                flint_printf("B:\n");
                fmpz_poly_mat_print(B, "x");
                flint_printf("X:\n");
                fmpz_poly_mat_print(X, "x");
                fflush(stdout);
                flint_abort();
            }
        }

        if (solved != !fmpz_poly_is_zero(den))
        {
            flint_printf("FAIL: return value does not match denominator\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_mul(AX, A, X);
        fmpz_poly_mat_scalar_mul_fmpz_poly(Bden, B, den);

        if (!fmpz_poly_mat_equal(AX, Bden))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            flint_printf("X:\n");
            fmpz_poly_mat_print(X, "x");
            flint_printf("AX:\n");
            fmpz_poly_mat_print(AX, "x");
            flint_printf("Bden:\n");
            fmpz_poly_mat_print(Bden, "x");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(den);
        fmpz_poly_clear(det);
        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(X);
        fmpz_poly_mat_clear(AX);
        fmpz_poly_mat_clear(Bden);
    }

    TEST_FUNCTION_END(state);
}
