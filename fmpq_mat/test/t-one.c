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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpq.h"
#include "fmpq_mat.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("one....");
    fflush(stdout);

    /* 1 * A == A * 1 == A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C, I;

        len_t n, bits;

        n = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, n);
        fmpq_mat_init(C, n, n);
        fmpq_mat_init(I, n, n);

        fmpq_mat_randtest(A, state, bits);
        fmpq_mat_one(I);

        fmpq_mat_mul(B, I, A);
        fmpq_mat_mul(C, A, I);

        result = fmpq_mat_equal(A, B) && fmpq_mat_equal(A, C);
        if (!result)
        {
            printf("FAIL:\n");
            printf("A:\n");
            fmpq_mat_print(A);
            printf("B:\n");
            fmpq_mat_print(B);
            printf("C:\n");
            fmpq_mat_print(C);
            printf("I:\n");
            fmpq_mat_print(I);
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        fmpq_mat_clear(I);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
