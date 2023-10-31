/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly_mat.h"

TEST_FUNCTION_START(nmod_poly_mat_solve_fflu, state)
{
    slong i;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, X, B, AX, Bden;
        nmod_poly_t den, det;
        slong n, m, deg;
        float density;
        int solved;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 15);
        m = n_randint(state, 5);
        deg = 1 + n_randint(state, 5);
        density = n_randint(state, 100) * 0.01;

        nmod_poly_mat_init(A, n, n, mod);
        nmod_poly_mat_init(B, n, m, mod);
        nmod_poly_mat_init(X, n, m, mod);
        nmod_poly_mat_init(AX, n, m, mod);
        nmod_poly_mat_init(Bden, n, m, mod);
        nmod_poly_init(den, mod);
        nmod_poly_init(det, mod);

        nmod_poly_mat_randtest_sparse(A, state, deg, density);
        nmod_poly_mat_randtest_sparse(B, state, deg, density);

        solved = nmod_poly_mat_solve_fflu(X, den, A, B);
        nmod_poly_mat_det_interpolate(det, A);

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
            if (!nmod_poly_equal(den, det))
            {
                nmod_poly_neg(det, det);
                flint_printf("FAIL: den != +/- det(A)\n");
                flint_printf("den:\n"); nmod_poly_print(den);
                flint_printf("\n\n");
                flint_printf("det:\n"); nmod_poly_print(det);
                flint_printf("\n\n");
                flint_printf("A:\n");
                nmod_poly_mat_print(A, "x");
                flint_printf("B:\n");
                nmod_poly_mat_print(B, "x");
                flint_printf("X:\n");
                nmod_poly_mat_print(X, "x");
                fflush(stdout);
                flint_abort();
            }
        }

        if (solved != !nmod_poly_is_zero(den))
        {
            flint_printf("FAIL: return value does not match denominator\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_mul(AX, A, X);
        nmod_poly_mat_scalar_mul_nmod_poly(Bden, B, den);

        if (!nmod_poly_mat_equal(AX, Bden))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            nmod_poly_mat_print(A, "x");
            flint_printf("B:\n");
            nmod_poly_mat_print(B, "x");
            flint_printf("X:\n");
            nmod_poly_mat_print(X, "x");
            flint_printf("AX:\n");
            nmod_poly_mat_print(AX, "x");
            flint_printf("Bden:\n");
            nmod_poly_mat_print(Bden, "x");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(den);
        nmod_poly_clear(det);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(X);
        nmod_poly_mat_clear(AX);
        nmod_poly_mat_clear(Bden);
    }

    TEST_FUNCTION_END(state);
}
