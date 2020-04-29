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
#include "flint.h"
#include "fmpz_sparse_mat.h"

int
main(void)
{
    slong rep, bits, r, c, i, j, k, nnz;
    fmpz_sparse_mat_t A, B, C;
    slong *rows;
    slong *cols;
    fmpz *vals;
    FLINT_TEST_INIT(state);
    
    flint_printf("construction from entries....");
    fflush(stdout);

    for (rep = 0; rep < 1000; rep++)
    {
        do bits = n_randint(state, 100);
        while (bits < UWORD(2));
        r = 1; /*n_randint(state, 10);*/
        c = n_randint(state, 10);
        fmpz_sparse_mat_init(A, r, c);
        fmpz_sparse_mat_init(B, r, c);
        fmpz_sparse_mat_init(C, 0, c);

        fmpz_sparse_mat_randtest(A, state, 0, c, bits);
        fmpz_sparse_mat_randtest(B, state, 0, c, bits);
        nnz = 0;
        for (i = 0; i < r; ++i) nnz += A->rows[i].nnz;

        /* Construct B from entries of A */
        rows = flint_malloc(nnz * sizeof(*rows));
        cols = flint_malloc(nnz * sizeof(*cols));
        vals = _fmpz_vec_init(nnz);
        for (i = k = 0; i < r; ++i)
        {
            for (j = 0; j < A->rows[i].nnz; ++j, ++k)
            {
                rows[k] = i;
                cols[k] = A->rows[i].entries[j].ind;
                fmpz_set(&vals[k], A->rows[i].entries[j].val);
            }
        }
        fmpz_sparse_mat_from_entries(B, rows, cols, vals, nnz);

        if (!fmpz_sparse_mat_equal(A, B))
        {
            flint_printf("FAIL: A != B\n");
            flint_printf("A = ");
            fmpz_sparse_mat_print_pretty(A);
            flint_printf("B = ");
            fmpz_sparse_mat_print_pretty(B);
            abort();
        }
        
        for (i = 0; i < r; ++i) fmpz_sparse_mat_append_row(C, &A->rows[i]);
        if (!fmpz_sparse_mat_equal(A, C))
        {
            flint_printf("FAIL: A != C\n");
            flint_printf("A = ");
            fmpz_sparse_mat_print_pretty(A);
            flint_printf("C = ");
            fmpz_sparse_mat_print_pretty(C);
            abort();
        }
        flint_free(rows);
        flint_free(cols);
        _fmpz_vec_clear(vals, nnz);
        fmpz_sparse_mat_clear(A);
        fmpz_sparse_mat_clear(B);
        fmpz_sparse_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
