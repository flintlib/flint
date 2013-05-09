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

    printf("sub....");
    fflush(stdout);

    flint_randinit(state);

    /* Check evaluation homomorphism */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        fmpz_mat_t a, b, c, d;
        fmpz_t x;
        len_t m, n, bits, deg;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, m, n);
        fmpz_poly_mat_init(C, m, n);

        fmpz_mat_init(a, m, n);
        fmpz_mat_init(b, m, n);
        fmpz_mat_init(c, m, n);
        fmpz_mat_init(d, m, n);

        fmpz_init(x);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);
        fmpz_poly_mat_sub(C, A, B);

        fmpz_randtest(x, state, 1 + n_randint(state, 100));

        fmpz_poly_mat_evaluate_fmpz(a, A, x);
        fmpz_poly_mat_evaluate_fmpz(b, B, x);
        fmpz_poly_mat_evaluate_fmpz(d, C, x);
        fmpz_mat_sub(c, a, b);

        if (!fmpz_mat_equal(c, d))
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

        fmpz_mat_clear(a);
        fmpz_mat_clear(b);
        fmpz_mat_clear(c);
        fmpz_mat_clear(d);

        fmpz_clear(x);
    }

    /* Check aliasing C and A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        len_t m, n, bits, deg;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, m, n);
        fmpz_poly_mat_init(C, m, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);

        fmpz_poly_mat_sub(C, A, B);
        fmpz_poly_mat_sub(A, A, B);

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
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, C;
        len_t m, n, bits, deg;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, m, n);
        fmpz_poly_mat_init(C, m, n);

        fmpz_poly_mat_randtest(A, state, deg, bits);
        fmpz_poly_mat_randtest(B, state, deg, bits);

        fmpz_poly_mat_sub(C, A, B);
        fmpz_poly_mat_sub(B, A, B);

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
