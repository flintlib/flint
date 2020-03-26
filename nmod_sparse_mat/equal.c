/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"

int
nmod_sparse_mat_equal(const nmod_sparse_mat_t mat1, const nmod_sparse_mat_t mat2)
{
    if (mat1->r != mat2->r) {
        printf("Different number of rows\n");
        return 0;
    }

    if (mat1->nnz != mat2->nnz) {
        printf("Different number of nonzeroes\n");
        return 0;
    }

    if (mat1->nnz == 0)
        return 1;

    if(memcmp(mat1->row_nnz, mat2->row_nnz, mat1->r * sizeof(*mat1->row_nnz))) { 
        printf("Different row_nnz\n"); return 0;
    }

    slong i, j;
    for(i=0; i<mat1->r; ++i) {
        nmod_sparse_mat_entry_struct *e1 = mat1->entries + mat1->row_starts[i];
        nmod_sparse_mat_entry_struct *e2 = mat2->entries + mat2->row_starts[i];
        for(j=0; j<mat1->row_nnz[i]; ++j, ++e1, ++e2) {
            if((e1->col - mat1->c_off != e2->col - mat2->c_off) || (e1->val != e2->val)) {
                flint_printf("Mismatched row %d, entry %d\n", i, j);
                return 0;
            }
        }
    }
    return 1;
}
