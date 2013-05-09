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

    printf("inv....");
    fflush(stdout);

    flint_randinit(state);

    /* Test aliasing */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, Ainv;
        nmod_poly_t den1, den2;
        len_t n, deg;
        float density;
        int ns1, ns2, result;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 8);
        deg = 1 + n_randint(state, 5);
        density = n_randint(state, 100) * 0.01;

        nmod_poly_mat_init(A, n, n, mod);
        nmod_poly_mat_init(Ainv, n, n, mod);
        nmod_poly_init(den1, mod);
        nmod_poly_init(den2, mod);

        nmod_poly_mat_randtest_sparse(A, state, deg, density);

        ns1 = nmod_poly_mat_inv(Ainv, den1, A);
        ns2 = nmod_poly_mat_inv(A, den2, A);

        result = ns1 == ns2;

        if (result && ns1 != 0)
        {
            result = nmod_poly_equal(den1, den2) &&
                nmod_poly_mat_equal(A, Ainv);
        }

        if (!result)
        {
            printf("FAIL (aliasing)!\n");
            nmod_poly_mat_print(A, "x"); printf("\n");
            nmod_poly_mat_print(Ainv, "x"); printf("\n");
            abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(Ainv);
        nmod_poly_clear(den1);
        nmod_poly_clear(den2);
    }

    /* Check A^(-1) = A = 1 */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, Ainv, B, Iden;
        nmod_poly_t den, det;
        len_t n, deg;
        float density;
        int nonsingular;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 10);
        deg = 1 + n_randint(state, 5);
        density = n_randint(state, 100) * 0.01;

        nmod_poly_mat_init(A, n, n, mod);
        nmod_poly_mat_init(Ainv, n, n, mod);
        nmod_poly_mat_init(B, n, n, mod);
        nmod_poly_mat_init(Iden, n, n, mod);
        nmod_poly_init(den, mod);
        nmod_poly_init(det, mod);

        nmod_poly_mat_randtest_sparse(A, state, deg, density);
        nonsingular = nmod_poly_mat_inv(Ainv, den, A);
        nmod_poly_mat_det_interpolate(det, A);

        if (n == 0)
        {
            if (nonsingular == 0 || !nmod_poly_is_one(den))
            {
                printf("FAIL: expected empty matrix to pass\n");
                abort();
            }
        }
        else
        {
            if (!nmod_poly_equal(den, det))
            {
                nmod_poly_neg(det, det);
                printf("FAIL: den != det(A)\n");
                abort();
            }

            nmod_poly_mat_mul(B, Ainv, A);
            nmod_poly_mat_one(Iden);
            nmod_poly_mat_scalar_mul_nmod_poly(Iden, Iden, den);

            if (!nmod_poly_mat_equal(B, Iden))
            {
                printf("FAIL:\n");
                printf("A:\n");
                nmod_poly_mat_print(A, "x");
                printf("Ainv:\n");
                nmod_poly_mat_print(Ainv, "x");
                printf("B:\n");
                nmod_poly_mat_print(B, "x");
                printf("den:\n");
                nmod_poly_print(den);
                abort();
            }
        }

        nmod_poly_clear(den);
        nmod_poly_clear(det);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(Ainv);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(Iden);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
