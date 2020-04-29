/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    slong rep, r, c, i, j, rk, *P, *Q;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_mat_t A, LU, L, U;
    nmod_mat_t dL, dU, dLU;
    FLINT_TEST_INIT(state);
    
    flint_printf("decomposing PAQ = LU....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        if (rep % 50 == 0) {flint_printf("."); fflush(stdout);}
        r = n_randint(state, 200);
        c = n_randint(state, 200);

        do n = n_randtest_not_zero(state);
        while (!n_is_prime(n));
        nmod_init(&mod, n);
        P = flint_malloc(r*sizeof(*P));
        Q = flint_malloc(c*sizeof(*P));
        nmod_sparse_mat_init(A, r, c, mod);
        nmod_sparse_mat_init(LU, r, c, mod);
        nmod_sparse_mat_init(L, r, c, mod);
        nmod_sparse_mat_init(U, r, c, mod);
        nmod_sparse_mat_randtest(A, state, 1, c);
        rk = nmod_sparse_mat_lu(P, Q, L, U, A);
        nmod_mat_init(dL, r, rk, n);
        nmod_mat_init(dU, rk, c, n);
        nmod_mat_init(dLU, r, c, n);
        
        /* Check that L is lower triangular (with ones on diagonal up to rank) */
        for (i = 0; i < r; ++i) 
        {
            if (i < rk && *nmod_sparse_vec_at(&L->rows[i], i) != UWORD(1))
            {
                flint_printf("FAIL: L does not have unit diagonal up to the rank\n");
            }
            for (j = 0; j < L->rows[i].nnz; ++j) 
            {
                nmod_sparse_entry_struct *e = &L->rows[i].entries[j];
                if (e->ind > i) 
                {
                    flint_printf("FAIL: L not lower triangular\n");
                    abort();
                }
                if (e->ind >= rk) 
                {
                    flint_printf("FAIL: L not trivial past the rank\n");
                    nmod_sparse_mat_print_pretty(L);
                    abort();
                }
            }
        }
        /* Check that U is upper triangular (with nonzero diagonal up to rank) */
        for (i = 0; i < r; ++i) 
        {
            if (i < rk && *nmod_sparse_vec_at(&U->rows[i], i) == UWORD(0))
            {
                flint_printf("FAIL: U does not have nonzero diagonal\n");
                abort();
            }
            if (i >= rk && U->rows[i].nnz != UWORD(0)) 
            {
                flint_printf("FAIL: U not trivial pas the rank\n");
                abort();
            }
            for (j = 0; j < U->rows[i].nnz; ++j) 
            {
                nmod_sparse_entry_struct *e = &U->rows[i].entries[j];
                if (e->ind < i) 
                {
                    flint_printf("FAIL: U not upper triangular\n");
                    abort();
                }
            }
        }
        nmod_sparse_mat_to_dense(dL, L);
        nmod_sparse_mat_to_dense(dU, U);
        nmod_mat_mul(dLU, dL, dU);
        nmod_sparse_mat_from_dense(LU, dLU);
        nmod_sparse_mat_permute_rows(A, P);
        nmod_sparse_mat_permute_cols(A, Q);
        if (!nmod_sparse_mat_equal(A, LU)) 
        {
            flint_printf("FAIL: PAQ != LU\n");
            flint_printf("PAQ=");
            nmod_sparse_mat_print_pretty(A);
            flint_printf("LU=");
            nmod_sparse_mat_print_pretty(LU);
            abort();
        }

        flint_free(P);
        flint_free(Q);
        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(U);
        nmod_sparse_mat_clear(L);
        nmod_sparse_mat_clear(LU);
        nmod_mat_clear(dL);
        nmod_mat_clear(dU);
        nmod_mat_clear(dLU);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
