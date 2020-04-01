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
    slong rep, r, c, i, j, k, nnz;
    mp_limb_t n;
    nmod_t mod;
    nmod_sparse_mat_t A, B, C;
    slong *rows;
    slong *cols;
    mp_limb_t *vals;
    FLINT_TEST_INIT(state);
    
    flint_printf("construction from entries....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        r = n_randint(state, 10);
        c = n_randint(state, 10);
        do n = n_randtest_not_zero(state);
        while (n == UWORD(1));
        nmod_init(&mod, n);
        nmod_sparse_mat_init(A, r, c, mod);
        nmod_sparse_mat_init(B, r, c, mod);
        nmod_sparse_mat_init(C, 0, c, mod);

        nmod_sparse_mat_randtest(A, state, 2, 2);
        nmod_sparse_mat_randtest(B, state, 2, 2);
        nnz = 0;
        for (i = 0; i < r; ++i) nnz += A->rows[i].nnz;

        /* Construct B from entries of A */
        rows = flint_malloc(nnz * sizeof(*rows));
        cols = flint_malloc(nnz * sizeof(*cols));
        vals = flint_malloc(nnz * sizeof(*vals));
        for (i = k = 0; i < r; ++i)
        {
            for (j = 0; j < A->rows[i].nnz; ++j, ++k)
            {
                rows[k] = i;
                cols[k] = A->rows[i].entries[j].ind;
                vals[k] = A->rows[i].entries[j].val;
            }
        }
        nmod_sparse_mat_from_entries(B, rows, cols, vals, nnz);

        if (!nmod_sparse_mat_equal(A, B))
        {
            flint_printf("FAIL: A != B\n");
            flint_printf("A = ");
            nmod_sparse_mat_print_pretty(A);
            flint_printf("B = ");
            nmod_sparse_mat_print_pretty(B);
            abort();
        }
        
        for (i = 0; i < r; ++i) nmod_sparse_mat_append_row(C, &A->rows[i]);
        if (!nmod_sparse_mat_equal(A, C))
        {
            flint_printf("FAIL: A != C\n");
            flint_printf("A = ");
            nmod_sparse_mat_print_pretty(A);
            flint_printf("C = ");
            nmod_sparse_mat_print_pretty(C);
            abort();
        }
        flint_free(rows);
        flint_free(cols);
        flint_free(vals);
        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(B);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
