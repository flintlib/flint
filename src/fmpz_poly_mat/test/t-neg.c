/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly_mat.h"

TEST_FUNCTION_START(fmpz_poly_mat_neg, state)
{
    slong i;

    /* Check evaluation homomorphism */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B;
        fmpz_mat_t a, b, c;
        fmpz_t x;
        slong m, n, bits, deg;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, m, n);

        fmpz_mat_init(a, m, n);
        fmpz_mat_init(b, m, n);
        fmpz_mat_init(c, m, n);

        fmpz_init(x);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_neg(B, A);

        fmpz_randtest(x, state, 1 + n_randint(state, 100));

        fmpz_poly_mat_evaluate_fmpz(a, A, x);
        fmpz_poly_mat_evaluate_fmpz(b, B, x);
        fmpz_mat_neg(c, a);

        if (!fmpz_mat_equal(b, c))
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

        fmpz_mat_clear(a);
        fmpz_mat_clear(b);
        fmpz_mat_clear(c);

        fmpz_clear(x);
    }

    /* Check aliasing B and A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B;
        slong m, n, bits, deg;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, m, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);

        fmpz_poly_mat_neg(B, A);
        fmpz_poly_mat_neg(A, A);

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
