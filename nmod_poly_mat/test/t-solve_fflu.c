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

    printf("solve....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, X, B, AX, Bden;
        nmod_poly_t den, det;
        len_t n, m, deg;
        float density;
        int solved;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        n = n_randint(state, 15);
        m = n_randint(state, 5);
        deg = 1 + n_randint(state, 5);
        density = n_randint(state, 100) * 0.01;

        nmod_poly_mat_init(A, n, n, mod);
        nmod_poly_mat_init(B, n, m, mod);
        nmod_poly_mat_init(X, n, m, mod);
        nmod_poly_mat_init(AX, n, m, mod);
        nmod_poly_mat_init(Bden, n, m, mod);
        nmod_poly_init(den, mod);
        nmod_poly_init(det, mod);

        nmod_poly_mat_randtest_sparse(A, state, deg, density);
        nmod_poly_mat_randtest_sparse(B, state, deg, density);

        solved = nmod_poly_mat_solve_fflu(X, den, A, B);
        nmod_poly_mat_det_interpolate(det, A);

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
            if (!nmod_poly_equal(den, det))
            {
                nmod_poly_neg(det, det);
                if (!nmod_poly_equal(den, det))
                {
                    nmod_poly_neg(det, det);
                    printf("FAIL: den != +/- det(A)\n");
                    printf("den:\n"); nmod_poly_print(den);
                    printf("\n\n");
                    printf("det:\n"); nmod_poly_print(det);
                    printf("\n\n");
                    printf("A:\n");
                    nmod_poly_mat_print(A, "x");
                    printf("B:\n");
                    nmod_poly_mat_print(B, "x");
                    printf("X:\n");
                    nmod_poly_mat_print(X, "x");
                    abort();
                }
            }
        }

        if (solved != !nmod_poly_is_zero(den))
        {
            printf("FAIL: return value does not match denominator\n");
            abort();
        }

        nmod_poly_mat_mul(AX, A, X);
        nmod_poly_mat_scalar_mul_nmod_poly(Bden, B, den);

        if (!nmod_poly_mat_equal(AX, Bden))
        {
            printf("FAIL:\n");
            printf("A:\n");
            nmod_poly_mat_print(A, "x");
            printf("B:\n");
            nmod_poly_mat_print(B, "x");
            printf("X:\n");
            nmod_poly_mat_print(X, "x");
            printf("AX:\n");
            nmod_poly_mat_print(AX, "x");
            printf("Bden:\n");
            nmod_poly_mat_print(Bden, "x");
            abort();
        }

        nmod_poly_clear(den);
        nmod_poly_clear(det);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(X);
        nmod_poly_mat_clear(AX);
        nmod_poly_mat_clear(Bden);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
