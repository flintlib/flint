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
#include "fmpz_mod_mat.h"
#include "ulong_extras.h"

void
fmpz_mod_mat_randrank(fmpz_mod_mat_t mat, flint_rand_t state, slong rank)
{
    slong i;
    fmpz * diag;

    if (rank < 0 || rank > mat->mat->r || rank > mat->mat->c)
    {
        flint_printf("Exception (fmpz_mod_mat_randrank). Impossible rank.\n");
        flint_abort();
    }

    diag = _fmpz_vec_init(rank);
    for (i = 0; i < rank; i++)
    {
        fmpz_randtest_mod(&diag[i], state, mat->mod);
        while (fmpz_is_zero(&diag[i]))
        {
            fmpz_randtest_mod(&diag[i], state, mat->mod);
        }
    }

    fmpz_mat_randpermdiag(mat->mat, state, diag, rank);

    _fmpz_vec_clear(diag, rank);
}

static void
check_rref(fmpz_mod_mat_t A)
{
    slong i, j, prev_pivot, prev_row_zero;

    prev_row_zero = 0;
    prev_pivot = -1;

    for (i = 0; i < A->mat->r; i++)
    {
        for (j = 0; j < A->mat->c; j++)
        {
            /* Found nonzero entry */
            if (!fmpz_is_zero(A->mat->rows[i] + j))
            {
                if (prev_row_zero)
                {
                    flint_printf("nonzero row after zero row\n");
                    abort();
                }

                if (j <= prev_pivot)
                {
                    flint_printf("pivot not strictly to the right of previous\n");
                    abort();
                }

                prev_pivot = j;
                break;
            }

            prev_row_zero = (j + 1 == A->mat->c);
        }
    }
}


int
main(void)
{
    fmpz_mod_mat_t A;
    fmpz_t p;
    slong i, m, n, d, r, rank;
    slong *perm;

    FLINT_TEST_INIT(state);

    flint_printf("rref....");
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
            fmpz_mod_mat_init(A, m, n, p);

            fmpz_mod_mat_randrank(A, state, r);

            rank = fmpz_mod_mat_rref(perm, A);

            if (r < rank)
            {
                fmpz_mod_mat_print_pretty(A);
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                abort();
            }

            check_rref(A);

            fmpz_mod_mat_clear(A);
        }

        fmpz_clear(p);
        flint_free(perm);
    }

    /* Dense */
    for (i = 0; i < 10000; i++)
    {
        m = n_randint(state, 5);
        n = n_randint(state, 4);
        perm = flint_malloc(FLINT_MAX(1, m) * sizeof(slong));

        fmpz_init_set_ui(p, n_randtest_prime(state, 0));

        for (r = 0; r <= FLINT_MIN(m, n); r++)
        {
            d = n_randint(state, 2 * m * n + 1);

            fmpz_mod_mat_init(A, m, n, p);
            fmpz_mod_mat_randrank(A, state, r);

            fmpz_mat_randops(A->mat, state, d);

            _fmpz_mod_mat_reduce(A);

            rank = fmpz_mod_mat_rref(perm, A);

            if (r < rank)
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong rank!\n");
                abort();
            }

            check_rref(A);

            fmpz_mod_mat_clear(A);
        }

        fmpz_clear(p);
        flint_free(perm);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
