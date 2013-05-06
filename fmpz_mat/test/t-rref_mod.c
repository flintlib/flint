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
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

static void
check_rref(fmpz_mat_t A)
{
    long i, j, prev_pivot, prev_row_zero;

    prev_row_zero = 0;
    prev_pivot = -1;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            /* Found nonzero entry */
            if (!fmpz_is_zero(A->rows[i] + j))
            {
                if (prev_row_zero)
                {
                    printf("nonzero row after zero row\n");
                    abort();
                }

                if (j <= prev_pivot)
                {
                    printf("pivot not strictly to the right of previous\n");
                    abort();
                }

                prev_pivot = j;
                break;
            }

            prev_row_zero = (j + 1 == A->c);
        }
    }
}


int
main(void)
{
    fmpz_mat_t A;
    fmpz_t p;
    flint_rand_t state;
    long i, j, k, m, n, b, d, r, rank;
    long *perm;

    printf("rref_mod....");
    fflush(stdout);

    flint_randinit(state);

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 10000; i++)
    {
        m = n_randint(state, 10);
        n = n_randint(state, 10);
        perm = flint_malloc(FLINT_MAX(1, m) * sizeof(long));

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);

            fmpz_mat_init(A, m, n);

            fmpz_mat_randrank(A, state, r, b);

            for (j = 0; j < m; j++)
                for (k = 0; k < n; k++)
                    fmpz_mod(fmpz_mat_entry(A, j, k), fmpz_mat_entry(A, j, k),
                             p);

            rank = fmpz_mat_rref_mod(perm, A, p);

            if (r < rank)
            {
                printf("FAIL:\n");
                printf("wrong rank!\n");
                abort();
            }

            check_rref(A);

            fmpz_mat_clear(A);
        }

        fmpz_clear(p);
        flint_free(perm);
    }

    /* Dense */
    for (i = 0; i < 10000; i++)
    {
        m = n_randint(state, 10);
        n = n_randint(state, 10);
        perm = flint_malloc(FLINT_MAX(1, m) * sizeof(long));

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2 * m * n + 1);

            fmpz_mat_init(A, m, n);
            fmpz_mat_randrank(A, state, r, b);

            fmpz_mat_randops(A, state, d);

            for (j = 0; j < m; j++)
                for (k = 0; k < n; k++)
                    fmpz_mod(fmpz_mat_entry(A, j, k), fmpz_mat_entry(A, j, k),
                             p);

            rank = fmpz_mat_rref_mod(perm, A, p);

            if (r < rank)
            {
                printf("FAIL:\n");
                printf("wrong rank!\n");
                abort();
            }

            check_rref(A);

            fmpz_mat_clear(A);
        }

        fmpz_clear(p);
        flint_free(perm);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
