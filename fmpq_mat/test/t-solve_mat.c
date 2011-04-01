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
#include <mpir.h>
#include "flint.h"
#include "fmpq.h"
#include "fmpq_mat.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("solve_mat....");
    fflush(stdout);

    for (i = 0; i < 1000; i++)
    {
        fmpq_mat_t A, B, X, AX;
        fmpq_t d;

        long n, m, bits;

        n = n_randint(state, 10);
        m = n_randint(state, 10);
        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, m);
        fmpq_mat_init(X, n, m);
        fmpq_mat_init(AX, n, m);

        fmpq_init(d);

        /* XXX: replace with a randtest function */
        do {
            fmpq_mat_randtest(A, state, bits);
            fmpq_mat_det(d, A);
        } while (fmpq_is_zero(d));

        fmpq_clear(d);

        fmpq_mat_solve_mat(X, A, B);
        fmpq_mat_mul(AX, A, X);

        if (!fmpq_mat_equal(AX, X))
        {
            printf("FAIL!\n");
            printf("A:\n");
            fmpq_mat_print(A);
            printf("B:\n");
            fmpq_mat_print(B);
            printf("X:\n");
            fmpq_mat_print(X);
            printf("AX:\n");
            fmpq_mat_print(AX);
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(X);
        fmpq_mat_clear(AX);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
