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
    len_t i;

    printf("det....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        fmpz_poly_t a, b, ab, c;
        len_t n, bits, deg;
        float density;

        n = n_randint(state, 10);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(B, n, n);
        fmpz_poly_mat_init(C, n, n);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(ab);
        fmpz_poly_init(c);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);
        fmpz_poly_mat_randtest_sparse(B, state, deg, bits, density);
        fmpz_poly_mat_mul(C, A, B);

        fmpz_poly_mat_det(a, A);
        fmpz_poly_mat_det(b, B);
        fmpz_poly_mat_det(c, C);
        fmpz_poly_mul(ab, a, b);

        if (!fmpz_poly_equal(c, ab))
        {
            printf("FAIL:\n");
            printf("determinants don't agree!\n");
            printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            printf("C:\n");
            fmpz_poly_mat_print(C, "x");
            printf("det(A):\n");
            fmpz_poly_print_pretty(a, "x");
            printf("\ndet(B):\n");
            fmpz_poly_print_pretty(b, "x");
            printf("\ndet(C):\n");
            fmpz_poly_print_pretty(c, "x");
            printf("\ndet(A)*det(B):\n");
            fmpz_poly_print_pretty(ab, "x");
            printf("\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(ab);
        fmpz_poly_clear(c);

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(C);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
