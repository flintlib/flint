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
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"


int
main(void)
{
    flint_rand_t state;
    long i;

    printf("sqr_KS....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, C, D;
        long m, bits, deg;

        /* TODO: add separate unsigned tests */
        m = n_randint(state, 15);
        deg = 1 + n_randint(state, 15);
        bits = 1 + n_randint(state, 150);

        fmpz_poly_mat_init(A, m, m);
        fmpz_poly_mat_init(C, m, m);
        fmpz_poly_mat_init(D, m, m);

        if (n_randint(state, 2))
            fmpz_poly_mat_randtest(A, state, deg, bits);
        else
            fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(C, state, deg, bits);  /* noise in output */

        fmpz_poly_mat_sqr_classical(C, A);
        fmpz_poly_mat_sqr_KS(D, A);

        if (!fmpz_poly_mat_equal(C, D))
        {
            printf("FAIL:\n");
            printf("products don't agree!\n");
            printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            printf("B:\n");
            fmpz_poly_mat_print(C, "x");
            printf("D:\n");
            fmpz_poly_mat_print(D, "x");
            printf("\n");
            abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(C);
        fmpz_poly_mat_clear(D);
    }

    /* Check aliasing B and A */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B;
        long m, bits, deg;

        m = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_poly_mat_init(A, m, m);
        fmpz_poly_mat_init(B, m, m);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);

        fmpz_poly_mat_sqr_KS(B, A);
        fmpz_poly_mat_sqr_KS(A, A);

        if (!fmpz_poly_mat_equal(B, A))
        {
            printf("FAIL (aliasing):\n");
            printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            printf("\n");
            abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
