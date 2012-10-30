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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "generics.h"
#include "perm.h"

int
main(void)
{
    flint_rand_t state;
    long i;

    printf("inv....");
    fflush(stdout);

    flint_randinit(state);

    /* Test aliasing */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        elem_mat_t A, Ainv;
        elem_poly_t den1, den2;
        ring_t ZZ, ZZx, MM;
        long n;
        long size[3] = {5, 5, 5};
        int ns1, ns2;
        int result;

        n = n_randint(state, 8);

        ring_init_fmpz(ZZ);
        ring_init_poly(ZZx, ZZ);
        ring_init_mat(MM, ZZx);

        elem_mat_init(A, n, n, MM);
        elem_mat_init(Ainv, n, n, MM);
        elem_init(den1, ZZx);
        elem_init(den2, ZZx);

        elem_mat_randtest(A, state, size, MM);

        ns1 = elem_mat_inv(Ainv, den1, A, MM);
        ns2 = elem_mat_inv(A, den2, A, MM);

        result = ns1 == ns2;

        if (result && ns1 != 0)
        {
            result = elem_equal(den1, den2, ZZx) &&
                elem_mat_equal(A, Ainv, MM);
        }

        if (!result)
        {
            printf("FAIL (aliasing)!\n");
            elem_mat_print(A, MM); printf("\n");
            elem_mat_print(Ainv, MM); printf("\n");
            abort();
        }

        elem_mat_clear(A, MM);
        elem_mat_clear(Ainv, MM);
        elem_clear(den1, ZZx);
        elem_clear(den2, ZZx);
    }

    /* Check A^(-1) = A = 1 */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        elem_mat_t A, Ainv, B, Iden;
        elem_poly_t den, det;
        ring_t ZZ, ZZx, MM;
        long n;
        long size[3] = {5, 5, 5};
        int nonsingular;

        n = n_randint(state, 8);

        ring_init_fmpz(ZZ);
        ring_init_poly(ZZx, ZZ);
        ring_init_mat(MM, ZZx);

        elem_mat_init(A, n, n, MM);
        elem_mat_init(Ainv, n, n, MM);
        elem_mat_init(B, n, n, MM);
        elem_mat_init(Iden, n, n, MM);
        elem_init(den, ZZx);
        elem_init(det, ZZx);

        elem_mat_randtest(A, state, size, MM);

        nonsingular = elem_mat_inv(Ainv, den, A, MM);
        elem_mat_det(det, A, MM);

        if (n == 0)
        {
            if (nonsingular == 0 || !elem_is_one(den, ZZx))
            {
                printf("FAIL: expected empty matrix to pass\n");
                abort();
            }
        }
        else
        {
/*
            if (!elem_equal(den, det, ZZx))
            {
                elem_neg(det, det, ZZx);
                printf("FAIL: den != det(A)\n");
                abort();
            }
*/

            if (nonsingular)
            {
                elem_mat_mul(B, Ainv, A, MM);
                elem_mat_one(Iden, MM);
                elem_mat_scalar_mul(Iden, Iden, den, MM);

                if (!elem_mat_equal(B, Iden, MM))
                {
                    printf("FAIL:\n");
                    printf("A:\n");
                    elem_mat_print(A, MM);
                    printf("Ainv:\n");
                    elem_mat_print(Ainv, MM);
                    printf("B:\n");
                    elem_mat_print(B, MM);
                    printf("den:\n");
                    elem_print(den, ZZx);
                    abort();
                }
            }
            else
            {
                if (!elem_is_zero(det, ZZx) || !elem_is_zero(den, ZZx))
                {
                    printf("FAIL (reported invertible):\n");
                    printf("A:\n");
                    elem_mat_print(A, MM);
                    printf("Ainv:\n");
                    elem_mat_print(Ainv, MM);
                    printf("den:\n");
                    elem_print(den, ZZx);
                    abort();
                }
            }
        }

        elem_clear(den, ZZx);
        elem_clear(det, ZZx);
        elem_mat_clear(A, MM);
        elem_mat_clear(Ainv, MM);
        elem_mat_clear(B, MM);
        elem_mat_clear(Iden, MM);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

