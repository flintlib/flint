/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

int
main(void)
{
    slong i;

    FLINT_TEST_INIT(state);

    flint_printf("mullow....");
    fflush(stdout);    

    /* Compare with mul */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C, D;
        slong m, n, k, bits, deg, len;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        k = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);
        len = n_randint(state, 10);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, n, k);
        fmpz_poly_mat_init(C, m, k);
        fmpz_poly_mat_init(D, m, k);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);
        fmpz_poly_mat_randtest(C, state, deg, bits);  /* noise in output */
        fmpz_poly_mat_randtest(D, state, deg, bits);  /* noise in output */

        fmpz_poly_mat_mullow(C, A, B, len);
        fmpz_poly_mat_mul(D, A, B);
        fmpz_poly_mat_truncate(D, len);

        if (!fmpz_poly_mat_equal(C, D))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            flint_printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            flint_printf("C:\n");
            fmpz_poly_mat_print(C, "x");
            flint_printf("D:\n");
            fmpz_poly_mat_print(D, "x");
            flint_printf("\n");
            abort();
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
        slong m, n, bits, deg, len;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);
        len = n_randint(state, 10);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, n, n);
        fmpz_poly_mat_init(C, m, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);
        fmpz_poly_mat_randtest(C, state, deg, bits);  /* noise in output */

        fmpz_poly_mat_mullow(C, A, B, len);
        fmpz_poly_mat_mullow(A, A, B, len);

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
            abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(C);
    }

    /* Check aliasing C and B */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        slong m, n, bits, deg, len;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);
        len = n_randint(state, 10);

        fmpz_poly_mat_init(A, m, m);
        fmpz_poly_mat_init(B, m, n);
        fmpz_poly_mat_init(C, m, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);
        fmpz_poly_mat_randtest(C, state, deg, bits);  /* noise in output */

        fmpz_poly_mat_mullow(C, A, B, len);
        fmpz_poly_mat_mullow(B, A, B, len);

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
            abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
