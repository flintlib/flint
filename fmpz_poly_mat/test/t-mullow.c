/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

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
    flint_rand_t state;
    len_t i;

    printf("mullow....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with mul */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C, D;
        len_t m, n, k, bits, deg, len;

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
            printf("FAIL:\n");
            printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            printf("C:\n");
            fmpz_poly_mat_print(C, "x");
            printf("D:\n");
            fmpz_poly_mat_print(D, "x");
            printf("\n");
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
        len_t m, n, bits, deg, len;

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
            printf("FAIL:\n");
            printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            printf("C:\n");
            fmpz_poly_mat_print(C, "x");
            printf("\n");
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
        len_t m, n, bits, deg, len;

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
            printf("FAIL:\n");
            printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            printf("C:\n");
            fmpz_poly_mat_print(C, "x");
            printf("\n");
            abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(C);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
