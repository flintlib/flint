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

TEST_FUNCTION_START(fmpz_poly_mat_pow_trunc, state)
{
    slong i;

    /* Compare with pow */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        slong n, exp, bits, deg, len;

        n = n_randint(state, 10);
        exp = n_randint(state, 15);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);
        len = n_randint(state, 10);

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(B, n, n);
        fmpz_poly_mat_init(C, n, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);  /* noise in output */
        fmpz_poly_mat_randtest(C, state, deg, bits);  /* noise in output */

        fmpz_poly_mat_pow_trunc(B, A, exp, len);
        fmpz_poly_mat_pow(C, A, exp);
        fmpz_poly_mat_truncate(C, len);

        if (!fmpz_poly_mat_equal(B, C))
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

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B;
        slong n, exp, bits, deg, len;

        n = n_randint(state, 10);
        exp = n_randint(state, 15);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);
        len = n_randint(state, 10);

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(B, n, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);

        fmpz_poly_mat_pow_trunc(B, A, exp, len);
        fmpz_poly_mat_pow_trunc(A, A, exp, len);

        if (!fmpz_poly_mat_equal(B, A))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
