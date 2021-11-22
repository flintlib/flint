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

    flint_printf("sqrlow....");
    fflush(stdout);    

    /* Compare with sqr */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        slong n, bits, deg, len;

        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);
        len = n_randint(state, 10);

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(B, n, n);
        fmpz_poly_mat_init(C, n, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);  /* noise in output */
        fmpz_poly_mat_randtest(C, state, deg, bits);  /* noise in output */

        fmpz_poly_mat_sqrlow(B, A, len);
        fmpz_poly_mat_sqr(C, A);
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
        slong n, bits, deg, len;

        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);
        len = n_randint(state, 10);

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(B, n, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);

        fmpz_poly_mat_sqrlow(B, A, len);
        fmpz_poly_mat_sqrlow(A, A, len);

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

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
