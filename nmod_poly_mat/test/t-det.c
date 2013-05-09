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
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "fmpz.h"

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
        nmod_poly_mat_t A, B, C;
        nmod_poly_t a, b, ab, c;
        len_t n, deg;
        mp_limb_t mod;
        float density;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 10);
        deg = 1 + n_randint(state, 5);
        density = n_randint(state, 100) * 0.01;

        nmod_poly_mat_init(A, n, n, mod);
        nmod_poly_mat_init(B, n, n, mod);
        nmod_poly_mat_init(C, n, n, mod);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);
        nmod_poly_init(ab, mod);
        nmod_poly_init(c, mod);

        nmod_poly_mat_randtest_sparse(A, state, deg, density);
        nmod_poly_mat_randtest_sparse(B, state, deg, density);
        nmod_poly_mat_mul(C, A, B);

        nmod_poly_mat_det(a, A);
        nmod_poly_mat_det(b, B);
        nmod_poly_mat_det(c, C);
        nmod_poly_mul(ab, a, b);

        if (!nmod_poly_equal(c, ab))
        {
            printf("FAIL:\n");
            printf("determinants don't agree!\n");
            printf("A:\n");
            nmod_poly_mat_print(A, "x");
            printf("B:\n");
            nmod_poly_mat_print(B, "x");
            printf("C:\n");
            nmod_poly_mat_print(C, "x");
            printf("det(A):\n");
            nmod_poly_print(a);
            printf("\ndet(B):\n");
            nmod_poly_print(b);
            printf("\ndet(C):\n");
            nmod_poly_print(c);
            printf("\ndet(A)*det(B):\n");
            nmod_poly_print(ab);
            printf("\n");
            abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(ab);
        nmod_poly_clear(c);

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
