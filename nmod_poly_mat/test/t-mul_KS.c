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

    printf("mul_KS....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, C, D;
        len_t m, n, k, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 15);
        n = n_randint(state, 15);
        k = n_randint(state, 15);
        deg = 1 + n_randint(state, 15);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(B, n, k, mod);
        nmod_poly_mat_init(C, m, k, mod);
        nmod_poly_mat_init(D, m, k, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_randtest(B, state, deg);
        nmod_poly_mat_randtest(C, state, deg);  /* noise in output */

        nmod_poly_mat_mul_classical(C, A, B);
        nmod_poly_mat_mul_KS(D, A, B);

        if (!nmod_poly_mat_equal(C, D))
        {
            printf("FAIL:\n");
            printf("products don't agree!\n");
            printf("A:\n");
            nmod_poly_mat_print(A, "x");
            printf("B:\n");
            nmod_poly_mat_print(B, "x");
            printf("C:\n");
            nmod_poly_mat_print(C, "x");
            printf("D:\n");
            nmod_poly_mat_print(D, "x");
            printf("\n");
            abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
        nmod_poly_mat_clear(D);
    }

    /* Check aliasing C and A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, C;
        len_t m, n, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(B, n, n, mod);
        nmod_poly_mat_init(C, m, n, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_randtest(B, state, deg);
        nmod_poly_mat_randtest(C, state, deg);  /* noise in output */

        nmod_poly_mat_mul_KS(C, A, B);
        nmod_poly_mat_mul_KS(A, A, B);

        if (!nmod_poly_mat_equal(C, A))
        {
            printf("FAIL:\n");
            printf("A:\n");
            nmod_poly_mat_print(A, "x");
            printf("B:\n");
            nmod_poly_mat_print(B, "x");
            printf("C:\n");
            nmod_poly_mat_print(C, "x");
            printf("\n");
            abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
    }

    /* Check aliasing C and B */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, C;
        len_t m, n, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        n = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(A, m, m, mod);
        nmod_poly_mat_init(B, m, n, mod);
        nmod_poly_mat_init(C, m, n, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_randtest(B, state, deg);
        nmod_poly_mat_randtest(C, state, deg);  /* noise in output */

        nmod_poly_mat_mul_KS(C, A, B);
        nmod_poly_mat_mul_KS(B, A, B);

        if (!nmod_poly_mat_equal(C, B))
        {
            printf("FAIL:\n");
            printf("A:\n");
            nmod_poly_mat_print(A, "x");
            printf("B:\n");
            nmod_poly_mat_print(B, "x");
            printf("C:\n");
            nmod_poly_mat_print(C, "x");
            printf("\n");
            abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(C);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
