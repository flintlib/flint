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

void
nmod_sparse_mat_window_init(nmod_sparse_mat_t B, const nmod_sparse_mat_t A, slong r1, slong c1, slong r2, slong c2)
{
    slong i, j;
    r1 = FLINT_MIN(r1, A->r);
    r2 = FLINT_MIN(r2, A->r);
    if(r2 < r1) r2 = r1;
    if(c2 < c1) c2 = c1;

    B->entries = A->entries;
    B->mod = A->mod;
    B->r = r2-r1;
    B->c_off = c1;
    B->c = c2-c1;
    B->nnz = 0;
    if(B->r == 0) {B->row_starts = B->row_nnz = NULL; return;}
    B->row_starts = flint_malloc(B->r * sizeof(*B->row_starts));
    B->row_nnz = flint_malloc(B->r * sizeof(*B->row_nnz));
    memcpy(B->row_starts, A->row_starts + r1, B->r*sizeof(*B->row_starts));
    memcpy(B->row_nnz, A->row_nnz + r1, B->r*sizeof(*B->row_nnz));

    for(i = r1; i<r2; ++i) {
        nmod_sparse_mat_entry_struct *Ae = A->entries + A->row_starts[i];
        for(j=0; j<A->row_nnz[i]; ++j) {
            if(c1 > 0 && Ae[j].col < c1) {
                B->row_starts[i-r1] += 1;
                B->row_nnz[i-r1] -= 1;
            } else if(Ae[j].col >= c2) {
                B->row_nnz[i-r1] -= 1;
            } else {
                B->nnz++;
            }
        }
    }
}
