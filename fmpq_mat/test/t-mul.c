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

    printf("mul....");
    fflush(stdout);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, C, D;

        long m, n, k, bits;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);

        bits = 1 + n_randint(state, 100);

        fmpq_mat_init(A, m, k);
        fmpq_mat_init(B, k, n);
        fmpq_mat_init(C, m, n);
        fmpq_mat_init(D, m, n);

        fmpq_mat_randtest(A, state, bits);
        fmpq_mat_randtest(B, state, bits);
        fmpq_mat_randtest(C, state, bits);  /* noise in output */

        fmpq_mat_mul_direct(C, A, B);
        fmpq_mat_mul_cleared(D, A, B);

        result = fmpq_mat_equal(C, D);
        if (!result)
        {
            printf("FAIL:\n");
            printf("A:\n");
            fmpq_mat_print(A);
            printf("B:\n");
            fmpq_mat_print(B);
            printf("C:\n");
            fmpq_mat_print(C);
            printf("D:\n");
            fmpq_mat_print(D);
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        fmpq_mat_clear(D);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
