/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
    slong i, j, prev_pivot, prev_row_zero;

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
                    flint_printf("nonzero row after zero row\n");
                    fflush(stdout);
                    flint_abort();
                }

                if (j <= prev_pivot)
                {
                    flint_printf("pivot not strictly to the right of previous\n");
                    fflush(stdout);
                    flint_abort();
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
    slong i, j, k, m, n, b, d, r, rank;
    slong *perm;

    FLINT_TEST_INIT(state);

    flint_printf("rref_mod....");
    fflush(stdout);    

    /* Maximally sparse matrices of given rank */
    for (i = 0; i < 10000; i++)
    {
        m = n_randint(state, 10);
        n = n_randint(state, 10);
        perm = flint_malloc(FLINT_MAX(1, m) * sizeof(slong));

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
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
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
        perm = flint_malloc(FLINT_MAX(1, m) * sizeof(slong));

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
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                fflush(stdout);
                flint_abort();
            }

            check_rref(A);

            fmpz_mat_clear(A);
        }

        fmpz_clear(p);
        flint_free(perm);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
