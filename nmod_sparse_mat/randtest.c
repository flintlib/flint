/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"

static int smat_ecmp(const void *va, const void *vb) {
    const nmod_sparse_mat_entry_struct *a = va;
    const nmod_sparse_mat_entry_struct *b = vb;
    if(a->col < b->col) return -1;
    if(b->col < a->col) return 1;
    return 0;
}
void
nmod_sparse_mat_randtest(nmod_sparse_mat_t mat, flint_rand_t state)
{
    slong i,j, entries_per_row = 0, m;

    if(mat->r >= 30) entries_per_row = 1+ n_randint(state, mat->r/10);
    if(entries_per_row < 2) entries_per_row = 2;
    mat->nnz = entries_per_row * mat->r;
    mat->entries = flint_malloc(mat->nnz*sizeof(*mat->entries));
    for (i = 0; i < mat->r; i++)
    {
        mat->row_starts[i] = i*entries_per_row;
        mat->row_nnz[i] = entries_per_row;
        nmod_sparse_mat_entry_struct *e = mat->entries + mat->row_starts[i];
        for(j = 0; j < entries_per_row; ++j) {
            e[j].col = j;
            do
            e[j].val = n_randtest(state) % mat->mod.n;
            while(e[j].val==UWORD(0));
        }

        /* Use resevoir sampling to get random support */
        for(j = entries_per_row; j < mat->r; ++j) {
            m = n_randint(state, j+1);
            if(m < entries_per_row) e[m].col = j;
        }
        qsort(e, entries_per_row, sizeof(*e), smat_ecmp);
    }
    _nmod_sparse_mat_set_c(mat);
}
