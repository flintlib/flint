/*
    Copyright (C) 2011 Fredrik Johansson

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

void nmod_sparse_mat_from_entries(nmod_sparse_mat_t mat, slong * rows, slong * cols, mp_limb_t * vals, slong nnz)
{
    slong r, i, j;
    for(r=i=0; r<mat->r; ++r, i=j) {
        mat->rows[r].nnz = 0;
        for(j=i; j<nnz && rows[j]==r; ++j);
        nmod_sparse_vec_from_entries(&mat->rows[r], cols+i, vals+i, j-i);
    }
}
