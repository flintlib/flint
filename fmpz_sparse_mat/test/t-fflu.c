/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "perm.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpz_sparse_mat.h"

int
main(void)
{
    slong rep, r, c, min_nnz, max_nnz, bits, rk, nreps = 1000, i, j;
    slong *P, *Q;
    fmpz_t *val;
    fmpz *D;
    fmpz_mat_t dL, dU, dLU;
    fmpz_sparse_mat_t L, U, M, LU;
    FLINT_TEST_INIT(state);

    flint_printf("fflu....");
    fflush(stdout);

    for (rep = 0; rep < nreps; rep++)
    {
        bits = 100; /*2 + n_randint(state, 10);*/
        r = c = n_randint(state, 10);
        /*c = n_randint(state, 20);*/
        min_nnz = 0;
        max_nnz = c;

        D = _fmpz_vec_init(r);
        P = flint_malloc(r*sizeof(*P));
        Q = flint_malloc(c*sizeof(*Q));
        fmpz_sparse_mat_init(M, r, c);
        fmpz_sparse_mat_init(L, r, c);
        fmpz_sparse_mat_init(U, r, c);
        fmpz_sparse_mat_init(LU, r, c);
        fmpz_sparse_mat_randtest(M, state, min_nnz, max_nnz, bits);
        rk = fmpz_sparse_mat_fflu(D, P, Q, L, U, M);

        /* Check that L is lower triangular (with ones on diagonal up to rank) */
        for (i = 0; i < r; ++i) 
        {
            val = fmpz_sparse_vec_at(&L->rows[i], i);
            if (i < rk && (val == NULL || !fmpz_is_one(*val)))
            {
                flint_printf("FAIL: L does not have unit diagonal up to the rank\n");
            }
            for (j = 0; j < L->rows[i].nnz; ++j) 
            {
                if (L->rows[i].entries[j].ind > i) 
                {
                    flint_printf("FAIL: L not lower triangular\n");
                    abort();
                }
                if (L->rows[i].entries[j].ind >= rk) 
                {
                    flint_printf("FAIL: L not trivial past the rank\n");
                    flint_printf("rank = %wd\n", rk);
                    flint_printf("L = ");
                    fmpz_sparse_mat_print_pretty(L);
                    abort();
                }
            }
        }
        /* Check that U is upper triangular (with nonzero diagonal up to rank) */
        for (i = 0; i < r; ++i) 
        {
            val = fmpz_sparse_vec_at(&U->rows[i], i);
            if (i < rk && (val == NULL || fmpz_is_zero(*val)))
            {
                flint_printf("FAIL: U does not have nonzero diagonal\n");
                abort();
            }
            if (i >= rk && U->rows[i].nnz != UWORD(0)) 
            {
                flint_printf("FAIL: U not trivial past the rank\n");
                abort();
            }
            for (j = 0; j < U->rows[i].nnz; ++j) 
            {
                if (U->rows[i].entries[j].ind < i) 
                {
                    flint_printf("FAIL: U not upper triangular\n");
                    abort();
                }
            }
        }

        
        fmpz_mat_init(dL, r, rk);
        fmpz_mat_init(dU, rk, c);
        fmpz_mat_init(dLU, r, c);

        /* Verify that PDMQ = LU */
        fmpz_sparse_mat_mul_diag_fmpz(M, M, D);
        fmpz_sparse_mat_permute_rows(M, P);
        fmpz_sparse_mat_permute_cols(M, Q);
        fmpz_sparse_mat_to_dense(dL, L);
        fmpz_sparse_mat_to_dense(dU, U);
        fmpz_mat_mul(dLU, dL, dU);
        fmpz_sparse_mat_from_dense(LU, dLU);
        if (!fmpz_sparse_mat_equal(M, LU)) 
        {
            flint_printf("FAIL: PDMQ != LU\n");
            flint_printf("PDMQ=");
            fmpz_sparse_mat_print_pretty(M);
            flint_printf("LU=");
            fmpz_sparse_mat_print_pretty(LU);
            flint_printf("diff = ");
            fmpz_sparse_mat_sub(LU, LU, M);
            fmpz_sparse_mat_print_pretty(LU);
            abort();
        }

        _fmpz_vec_clear(D, r);
        flint_free(P);
        flint_free(Q);
        fmpz_sparse_mat_clear(M);
        fmpz_sparse_mat_clear(U);
        fmpz_sparse_mat_clear(L);
        fmpz_sparse_mat_clear(LU);
        fmpz_mat_clear(dL);
        fmpz_mat_clear(dU);
        fmpz_mat_clear(dLU);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

