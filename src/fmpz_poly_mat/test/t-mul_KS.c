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

TEST_FUNCTION_START(fmpz_poly_mat_mul_KS, state)
{
    slong i;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C, D;
        slong m, n, k, bits, deg;

        m = n_randint(state, 15);
        n = n_randint(state, 15);
        k = n_randint(state, 15);
        deg = 1 + n_randint(state, 15);
        bits = 1 + n_randint(state, 150);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, n, k);
        fmpz_poly_mat_init(C, m, k);
        fmpz_poly_mat_init(D, m, k);

        if (n_randint(state, 2))
            fmpz_poly_mat_randtest(A, state, deg, bits);
        else
            fmpz_poly_mat_randtest_unsigned(A, state, deg, bits);

        if (n_randint(state, 2))
            fmpz_poly_mat_randtest(B, state, deg, bits);
        else
            fmpz_poly_mat_randtest_unsigned(B, state, deg, bits);

        fmpz_poly_mat_randtest(C, state, deg, bits);  /* noise in output */

        fmpz_poly_mat_mul_classical(C, A, B);
        fmpz_poly_mat_mul_KS(D, A, B);

        if (!fmpz_poly_mat_equal(C, D))
        {
            flint_printf("FAIL:\n");
            flint_printf("products don't agree!\n");
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            flint_printf("C:\n");
            fmpz_poly_mat_print(C, "x");
            flint_printf("D:\n");
            fmpz_poly_mat_print(D, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(C);
        fmpz_poly_mat_clear(D);
    }

    /* Check aliasing C and A */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        slong m, n, bits, deg;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, n, n);
        fmpz_poly_mat_init(C, m, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);
        fmpz_poly_mat_randtest(C, state, deg, bits);  /* noise in output */

        fmpz_poly_mat_mul_KS(C, A, B);
        fmpz_poly_mat_mul_KS(A, A, B);

        if (!fmpz_poly_mat_equal(C, A))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            flint_printf("C:\n");
            fmpz_poly_mat_print(C, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(C);
    }

    /* Check aliasing C and B */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        slong m, n, bits, deg;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_poly_mat_init(A, m, m);
        fmpz_poly_mat_init(B, m, n);
        fmpz_poly_mat_init(C, m, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);
        fmpz_poly_mat_randtest(C, state, deg, bits);  /* noise in output */

        fmpz_poly_mat_mul_KS(C, A, B);
        fmpz_poly_mat_mul_KS(B, A, B);

        if (!fmpz_poly_mat_equal(C, B))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            flint_printf("C:\n");
            fmpz_poly_mat_print(C, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
