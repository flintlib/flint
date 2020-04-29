/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"
#include "hashmap.h"
#include "longlong.h"

void _fmpz_sparse_mat_with_transpose_fix_support(fmpz_sparse_mat_with_transpose_t MT, slong r, slong *osupp, slong onnz)
{
    slong i = 0, j = 0, oc, nc, c;
    fmpz_sparse_vec_struct *row = &MT->M->rows[r];

    i = 0;
    while (1)
    {
        oc = (i==onnz) ? MT->M->c : osupp[i];
        nc = (j==row->nnz) ? MT->M->c : row->entries[j].ind;
        c = FLINT_MIN(oc, nc);
        if (c >= MT->M->c) break;
        if (oc < nc) hashmap_rem(&MT->cols[c], r);
        else hashmap_put(&MT->cols[c], r, &MT->M->rows[r].entries[j].val);
        if (oc <= nc) ++i;
        if (oc >= nc) ++j;
    }
}

