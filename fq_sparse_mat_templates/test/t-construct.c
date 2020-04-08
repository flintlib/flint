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
    slong rep, r, c, i, j, k, nnz;
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, sparse_mat_t) A, B, C;
    slong *rows;
    slong *cols;
    TEMPLATE(T, struct) *vals;
    FLINT_TEST_INIT(state);
    
    flint_printf("construction from entries....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        TEMPLATE(T, ctx_randtest) (ctx, state);
        r = n_randint(state, 10);
        c = n_randint(state, 10);
        TEMPLATE(T, sparse_mat_init) (A, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (B, r, c, ctx);
        TEMPLATE(T, sparse_mat_init) (C, 0, c, ctx);

        TEMPLATE(T, sparse_mat_randtest) (A, state, 2, 2, ctx);
        TEMPLATE(T, sparse_mat_randtest) (B, state, 2, 2, ctx);
        nnz = 0;
        for (i = 0; i < r; ++i) nnz += A->rows[i].nnz;

        /* Construct B from entries of A */
        rows = flint_malloc(nnz * sizeof(*rows));
        cols = flint_malloc(nnz * sizeof(*cols));
        vals = _TEMPLATE(T, vec_init) (nnz, ctx);
        for (i = k = 0; i < r; ++i)
        {
            for (j = 0; j < A->rows[i].nnz; ++j, ++k)
            {
                rows[k] = i;
                cols[k] = A->rows[i].entries[j].ind;
                TEMPLATE(T, set) (&vals[k], A->rows[i].entries[j].val, ctx);
            }
        }
        TEMPLATE(T, sparse_mat_from_entries) (B, rows, cols, vals, nnz, ctx);

        if (!TEMPLATE(T, sparse_mat_equal) (A, B, ctx))
        {
            flint_printf("FAIL: A != B\n");
            flint_printf("A = ");
            TEMPLATE(T, sparse_mat_print_pretty) (A, ctx);
            flint_printf("B = ");
            TEMPLATE(T, sparse_mat_print_pretty) (B, ctx);
            abort();
        }
        
        for (i = 0; i < r; ++i) TEMPLATE(T, sparse_mat_append_row) (C, &A->rows[i], ctx);
        if (!TEMPLATE(T, sparse_mat_equal) (A, C, ctx))
        {
            flint_printf("FAIL: A != C\n");
            flint_printf("A = ");
            TEMPLATE(T, sparse_mat_print_pretty) (A, ctx);
            flint_printf("C = ");
            TEMPLATE(T, sparse_mat_print_pretty) (C, ctx);
            abort();
        }
        flint_free(rows);
        flint_free(cols);
        _TEMPLATE(T, vec_clear) (vals, nnz, ctx);
        TEMPLATE(T, sparse_mat_clear) (A, ctx);
        TEMPLATE(T, sparse_mat_clear) (B, ctx);
        TEMPLATE(T, sparse_mat_clear) (C, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

#endif
