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
    slong m, mod, rep;
    FLINT_TEST_INIT(state);
    
    flint_printf("conversion to/from dense matrix....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        slong i, j, k;
        nmod_sparse_mat_t A, B, C;
        slong *rows, *cols;
        mp_limb_t *vals;
        nmod_sparse_mat_entry_struct *Ae;

        m = n_randint(state, 200);

        do
            mod = n_randtest_not_zero(state);
        while(mod <= 1);

        nmod_sparse_mat_init(A, m, mod);
        nmod_sparse_mat_init(B, m, mod);
        nmod_sparse_mat_init(C, m, mod);
        nmod_sparse_mat_randtest(A, state);
        
        /* Construct B from entries of A */
        rows = flint_malloc(A->nnz * sizeof(*rows));
        cols = flint_malloc(A->nnz * sizeof(*cols));
        vals = flint_malloc(A->nnz * sizeof(*vals));
        for(i=0, k=0; i<A->r; ++i) {
            for(j=0; j<A->row_nnz[i]; ++j, ++k) {
                rows[k] = i;
                cols[k] = A->entries[k].col;
                vals[k] = A->entries[k].val;
            }
        }
        nmod_sparse_mat_set_from_entries(B, rows, cols, vals, A->nnz);

        if (!nmod_sparse_mat_equal(A, B))
        {
            flint_printf("FAIL: A != B\n");
            abort();
        }
        /* Construct C from rows of A */
        for(i=0; i<A->r; ++i) {
            nmod_sparse_mat_append_row(C, i, cols + A->row_starts[i], vals + A->row_starts[i], A->row_nnz[i]);
        }

        if (!nmod_sparse_mat_equal(A, C))
        {
            flint_printf("FAIL: A != C\n");
            abort();
        }
        flint_free(rows);
        flint_free(cols);
        flint_free(vals);
        nmod_sparse_mat_clear(A);
        nmod_sparse_mat_clear(B);
        nmod_sparse_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
