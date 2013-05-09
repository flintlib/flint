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

    printf("sqrlow....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with sqr */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        len_t n, bits, deg, len;

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

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B;
        len_t n, bits, deg, len;

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
            printf("FAIL:\n");
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
