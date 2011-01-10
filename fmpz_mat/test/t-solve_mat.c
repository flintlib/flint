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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    fmpz_mat_t A, X, B, AX;
    fmpz_t den;
    flint_rand_t rnd;
    long i, m, n, r;

    printf("solve_mat....");
    fflush(stdout);

    fmpz_randinit(rnd);

    for (i = 0; i < 10000; i++)
    {
        m = n_randint(10);
        n = n_randint(10);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, m, n);
        fmpz_mat_init(AX, m, n);
        fmpz_init(den);

        fmpz_mat_randrank(A, rnd, m, 1+n_randint(2)*n_randint(100));
        fmpz_mat_randtest(B, rnd, 1+n_randint(2)*n_randint(100));

        /* Dense */
        if (n_randint(2))
            fmpz_mat_randops(A, rnd, 1+n_randint(1+m*m));

        fmpz_mat_solve_mat(X, den, A, B);

        fmpz_mat_mul(AX, A, X);
        _fmpz_vec_scalar_divexact_fmpz(AX->entries, AX->entries, m*n, den);

        if (!fmpz_mat_equal(AX, B))
        {
            printf("FAIL:\n");
            printf("AX != B!\n");
            printf("A:\n"),      fmpz_mat_print_pretty(A),  printf("\n");
            printf("B:\n"),      fmpz_mat_print_pretty(B),  printf("\n");
            printf("X:\n"),      fmpz_mat_print_pretty(X),  printf("\n");
            printf("den(X) = "), fmpz_print(den),           printf("\n");
            printf("AX:\n"),     fmpz_mat_print_pretty(AX), printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_mat_clear(AX);
        fmpz_clear(den);
    }

    /* Test singular systems */
    for (i = 0; i < 10000; i++)
    {
        m = 1 + n_randint(10);
        n = 1 + n_randint(10);
        r = n_randint(m);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, n);
        fmpz_mat_init(X, m, n);
        fmpz_mat_init(AX, m, n);
        fmpz_init(den);

        fmpz_mat_randrank(A, rnd, r, 1+n_randint(2)*n_randint(100));
        fmpz_mat_randtest(B, rnd, 1+n_randint(2)*n_randint(100));

        /* Dense */
        if (n_randint(2))
            fmpz_mat_randops(A, rnd, 1+n_randint(1+m*m));

        fmpz_mat_solve_mat(X, den, A, B);

        if (!fmpz_is_zero(den))
        {
            printf("FAIL:\n");
            printf("singular system gave nonzero determinant\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(X);
        fmpz_mat_clear(AX);
        fmpz_clear(den);
    }

    fmpz_mat_randclear(rnd);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
