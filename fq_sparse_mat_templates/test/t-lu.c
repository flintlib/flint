/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>
#include "ulong_extras.h"

int
main(void)
{
    slong rep, r, c, i, j, rk, *P, *Q;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, t) *val;
    TEMPLATE(T, sparse_mat_t) A, LU, L, U;
    TEMPLATE(T, mat_t) dL, dU, dLU;
    FLINT_TEST_INIT(state);
    
    flint_printf("decomposing PAQ = LU....");
    fflush(stdout);

    for (rep = 0; rep < 200; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        if (rep % 20 == 0) {flint_printf("."); fflush(stdout);}
        r = n_randint(state, 100);
        c = n_randint(state, 100);

        P = flint_malloc(r*sizeof(*P));
        Q = flint_malloc(c*sizeof(*P));
        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (LU, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (L, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (U, r, c, ctx);
        TEMPLATE(T, sparse_mat_randtest) (A, state, 1, c, ctx);
        rk = TEMPLATE(T, sparse_mat_lu) (P, Q, L, U, A, ctx);
        TEMPLATE(T, mat_init) (dL, r, rk, ctx);
        TEMPLATE(T, mat_init) (dU, rk, c, ctx);
        TEMPLATE(T, mat_init) (dLU, r, c, ctx);
        
        /* Check that L is lower triangular (with ones on diagonal up to rank) */
        for (i = 0; i < r; ++i) 
        {
            val = TEMPLATE(T, sparse_vec_at) (&L->rows[i], i, ctx);
            if (i < rk && (val == NULL || !TEMPLATE(T, is_one) (*val, ctx)))
            {
                flint_printf("FAIL: L does not have unit diagonal up to the rank\n");
            }
            for (j = 0; j < L->rows[i].nnz; ++j) 
            {
                TEMPLATE(T, sparse_entry_struct) *e = &L->rows[i].entries[j];
                if (e->ind > i) 
                {
                    flint_printf("FAIL: L not lower triangular\n");
                    abort();
                }
                if (e->ind >= rk) 
                {
                    flint_printf("FAIL: L not trivial past the rank\n");
                    /*TEMPLATE(T, sparse_mat_print_pretty) (L, ctx);*/
                    abort();
                }
            }
        }
        /* Check that U is upper triangular (with nonzero diagonal up to rank) */
        for (i = 0; i < r; ++i) 
        {
            val = TEMPLATE(T, sparse_vec_at) (&U->rows[i], i, ctx);
            if (i < rk && (val == NULL || TEMPLATE(T, is_zero) (*val, ctx)))
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
                TEMPLATE(T, sparse_entry_struct) *e = &U->rows[i].entries[j];
                if (e->ind < i) 
                {
                    flint_printf("FAIL: U not upper triangular\n");
                    abort();
                }
            }
        }
        TEMPLATE(T, sparse_mat_to_dense) (dL, L, ctx);
        TEMPLATE(T, sparse_mat_to_dense) (dU, U, ctx);
        TEMPLATE(T, mat_mul) (dLU, dL, dU, ctx);
        TEMPLATE(T, sparse_mat_from_dense) (LU, dLU, ctx);
        TEMPLATE(T, sparse_mat_permute_rows) (A, P, ctx);
        TEMPLATE(T, sparse_mat_permute_cols) (A, Q, ctx);
        if (!TEMPLATE(T, sparse_mat_equal) (A, LU, ctx)) 
        {
            flint_printf("FAIL: PAQ != LU\n");
            flint_printf("PAQ=");
            TEMPLATE(T, sparse_mat_print_pretty) (A, ctx);
            flint_printf("LU=");
            TEMPLATE(T, sparse_mat_print_pretty) (LU, ctx);
            abort();
        }

        flint_free(P);
        flint_free(Q);
        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, sparse_mat_clear) (U, ctx);
        TEMPLATE(T, sparse_mat_clear) (L, ctx);
        TEMPLATE(T, sparse_mat_clear) (LU, ctx);
        TEMPLATE(T, mat_clear) (dL, ctx);
        TEMPLATE(T, mat_clear) (dU, ctx);
        TEMPLATE(T, mat_clear) (dLU, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif
