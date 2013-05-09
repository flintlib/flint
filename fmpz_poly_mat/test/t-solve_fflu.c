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

    printf("solve....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, X, B, AX, Bden;
        fmpz_poly_t den, det;
        len_t n, m, bits, deg;
        float density;
        int solved;

        n = n_randint(state, 15);
        m = n_randint(state, 5);
        deg = 1 + n_randint(state, 5);
        bits = 1 + n_randint(state, 100);
        density = n_randint(state, 100) * 0.01;

        fmpz_poly_mat_init(A, n, n);
        fmpz_poly_mat_init(B, n, m);
        fmpz_poly_mat_init(X, n, m);
        fmpz_poly_mat_init(AX, n, m);
        fmpz_poly_mat_init(Bden, n, m);
        fmpz_poly_init(den);
        fmpz_poly_init(det);

        fmpz_poly_mat_randtest_sparse(A, state, deg, bits, density);
        fmpz_poly_mat_randtest_sparse(B, state, deg, bits, density);

        solved = fmpz_poly_mat_solve_fflu(X, den, A, B);
        fmpz_poly_mat_det_interpolate(det, A);

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
            if (!fmpz_poly_equal(den, det))
            {
                fmpz_poly_neg(det, det);
                if (!fmpz_poly_equal(den, det))
                {
                    fmpz_poly_neg(det, det);
                    printf("FAIL: den != +/- det(A)\n");
                    printf("den:\n"); fmpz_poly_print_pretty(den, "x");
                    printf("\n\n");
                    printf("det:\n"); fmpz_poly_print_pretty(det, "x");
                    printf("\n\n");
                    printf("A:\n");
                    fmpz_poly_mat_print(A, "x");
                    printf("B:\n");
                    fmpz_poly_mat_print(B, "x");
                    printf("X:\n");
                    fmpz_poly_mat_print(X, "x");
                    abort();
                }
            }
        }

        if (solved != !fmpz_poly_is_zero(den))
        {
            printf("FAIL: return value does not match denominator\n");
            abort();
        }

        fmpz_poly_mat_mul(AX, A, X);
        fmpz_poly_mat_scalar_mul_fmpz_poly(Bden, B, den);

        if (!fmpz_poly_mat_equal(AX, Bden))
        {
            printf("FAIL:\n");
            printf("A:\n");
            fmpz_poly_mat_print(A, "x");
            printf("B:\n");
            fmpz_poly_mat_print(B, "x");
            printf("X:\n");
            fmpz_poly_mat_print(X, "x");
            printf("AX:\n");
            fmpz_poly_mat_print(AX, "x");
            printf("Bden:\n");
            fmpz_poly_mat_print(Bden, "x");
            abort();
        }

        fmpz_poly_clear(den);
        fmpz_poly_clear(det);
        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(X);
        fmpz_poly_mat_clear(AX);
        fmpz_poly_mat_clear(Bden);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
