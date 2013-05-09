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
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("rref....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        len_t m, n, r, rank, b, d;
        fmpq_mat_t A, B, C;
        fmpz_mat_t M;
        fmpz_t den;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_init(den);

        for (r = 0; r <= FLINT_MIN(m,n); r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2*m*n + 1);

            fmpz_mat_init(M, m, n);
            fmpq_mat_init(A, m, n);
            fmpq_mat_init(B, m, n);
            fmpq_mat_init(C, m, n);

            fmpz_mat_randrank(M, state, r, b);

            if (i % 2 == 0)
                fmpz_mat_randops(M, state, d);

            fmpz_randtest_not_zero(den, state, b);
            fmpq_mat_set_fmpz_mat_div_fmpz(A, M, den);

            rank = fmpq_mat_rref_classical(B, A);
            if (r != rank)
            {
                printf("FAIL:\n");
                printf("fmpq_mat_rref_classical: wrong rank!\n");
                fmpq_mat_print(A);
                printf("\nrank: %ld computed: %ld\n", r, rank);
                abort();
            }

            rank = fmpq_mat_rref_fraction_free(C, A);
            if (r != rank)
            {
                printf("FAIL:\n");
                printf("fmpq_mat_rref_fraction_free: wrong rank!\n");
                abort();
            }

            if (!fmpq_mat_equal(B, C))
            {
                printf("FAIL:\n");
                printf("different results!\n");
                printf("A:\n");
                fmpq_mat_print(A);
                printf("\nB:\n");
                fmpq_mat_print(B);
                printf("\nC:\n");
                fmpq_mat_print(C);
                abort();
            }

            fmpz_mat_clear(M);
            fmpq_mat_clear(A);
            fmpq_mat_clear(B);
            fmpq_mat_clear(C);
        }

        fmpz_clear(den);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
