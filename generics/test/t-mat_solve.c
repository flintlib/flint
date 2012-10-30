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
    long iter;
    flint_rand_t state;

    printf("mat_solve....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        elem_mat_t A, X, B, AX, Bden;
        elem_poly_t den, det;
        ring_t ZZ, ZZx, MM;
        long m, n;
        int solved;
        long size[3] = {5, 5, 5};

        m = n_randint(state, 7);
        n = n_randint(state, 7);

        ring_init_fmpz(ZZ);
        ring_init_poly(ZZx, ZZ);
        ring_init_mat(MM, ZZx);

        elem_mat_init(A, n, n, MM);
        elem_mat_init(B, n, m, MM);
        elem_mat_init(X, n, m, MM);
        elem_mat_init(AX, n, m, MM);
        elem_mat_init(Bden, n, m, MM);
        elem_init(den, ZZx);
        elem_init(det, ZZx);

        elem_mat_randtest(A, state, size, MM);
        elem_mat_randtest(B, state, size, MM);

        solved = elem_mat_solve(X, den, A, B, MM);
        elem_mat_det(det, A, MM);

        if (m == 0 || n == 0)
        {
            if (solved == 0)
            {
                printf("FAIL: expected empty system to pass\n");
                abort();
            }
        }
        else
        {
            if (!elem_equal(den, det, ZZx))
            {
                elem_neg(det, det, ZZx);
                if (!elem_equal(den, det, ZZx))
                {
                    elem_neg(det, det, ZZx);
                    printf("FAIL: den != +/- det(A)\n");
                    printf("den:\n"); elem_print(den, ZZx);
                    printf("\n\n");
                    printf("det:\n"); elem_print(det, ZZx);
                    printf("\n\n");
                    printf("A:\n");
                    elem_mat_print(A, MM);
                    printf("B:\n");
                    elem_mat_print(B, MM);
                    printf("X:\n");
                    elem_mat_print(X, MM);
                    abort();
                }
            }
        }

        if (solved != !elem_is_zero(den, ZZx))
        {
            printf("FAIL: return value does not match denominator\n");
            abort();
        }

        elem_mat_mul(AX, A, X, MM);
        elem_mat_scalar_mul(Bden, B, den, MM);

        if (!elem_mat_equal(AX, Bden, MM))
        {
            printf("FAIL:\n");
            printf("A:\n");
            elem_mat_print(A, MM);
            printf("B:\n");
            elem_mat_print(B, MM);
            printf("X:\n");
            elem_mat_print(X, MM);
            printf("AX:\n");
            elem_mat_print(AX, MM);
            printf("Bden:\n");
            elem_mat_print(Bden, MM);
            abort();
        }

        elem_clear(den, ZZx);
        elem_clear(det, ZZx);
        elem_mat_clear(A, MM);
        elem_mat_clear(B, MM);
        elem_mat_clear(X, MM);
        elem_mat_clear(AX, MM);
        elem_mat_clear(Bden, MM);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

