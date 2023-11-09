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

TEST_FUNCTION_START(fmpz_poly_mat_inv, state)
{
    slong i;

    /* Test aliasing */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, Ainv;
        fmpz_poly_t den1, den2;
        slong n, bits, deg;
        float density;
        int ns1, ns2;
        int result;

        n = n_randint(state, 8);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(Ainv, n, n);
        fmpz_poly_init(den1);
        fmpz_poly_init(den2);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);

        ns1 = fmpz_poly_mat_inv(Ainv, den1, A);
        ns2 = fmpz_poly_mat_inv(A, den2, A);

        result = ns1 == ns2;

        if (result && ns1 != 0)
        {
            result = fmpz_poly_equal(den1, den2) &&
                fmpz_poly_mat_equal(A, Ainv);
        }

        if (!result)
        {
            flint_printf("FAIL (aliasing)!\n");
            fmpz_poly_mat_print(A, "x"); flint_printf("\n");
            fmpz_poly_mat_print(Ainv, "x"); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(Ainv);
        fmpz_poly_clear(den1);
        fmpz_poly_clear(den2);
    }

    /* Check A^(-1) = A = 1 */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, Ainv, B, Iden;
        fmpz_poly_t den, det;
        slong n, bits, deg;
        float density;
        int nonsingular;

        n = n_randint(state, 10);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(Ainv, n, n);
        fmpz_poly_mat_init(B, n, n);
        fmpz_poly_mat_init(Iden, n, n);
        fmpz_poly_init(den);
        fmpz_poly_init(det);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);
        nonsingular = fmpz_poly_mat_inv(Ainv, den, A);
        fmpz_poly_mat_det_interpolate(det, A);

        if (n == 0)
        {
            if (nonsingular == 0 || !fmpz_poly_is_one(den))
            {
                flint_printf("FAIL: expected empty matrix to pass\n");
                fflush(stdout);
                flint_abort();
            }
        }
        else
        {
            if (!fmpz_poly_equal(den, det))
            {
                fmpz_poly_neg(det, det);
                flint_printf("FAIL: den != det(A)\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_poly_mat_mul(B, Ainv, A);
            fmpz_poly_mat_one(Iden);
            fmpz_poly_mat_scalar_mul_fmpz_poly(Iden, Iden, den);

            if (!fmpz_poly_mat_equal(B, Iden))
            {
                flint_printf("FAIL:\n");
                flint_printf("A:\n");
                fmpz_poly_mat_print(A, "x");
                flint_printf("Ainv:\n");
                fmpz_poly_mat_print(Ainv, "x");
                flint_printf("B:\n");
                fmpz_poly_mat_print(B, "x");
                flint_printf("den:\n");
                fmpz_poly_print_pretty(den, "x");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_poly_clear(den);
        fmpz_poly_clear(det);
        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(Ainv);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(Iden);
    }

    TEST_FUNCTION_END(state);
}
