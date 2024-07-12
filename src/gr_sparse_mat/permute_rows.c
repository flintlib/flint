/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

int gr_lil_mat_permute_rows(gr_lil_mat_t mat, const slong * perm, gr_ctx_t ctx)
{
    slong i, j;
    slong *iperm;

    /* todo: bounds checking */
    if (perm == NULL)
    {
        return GR_DOMAIN;
    }

    // Get inverse permutation, i.e., iperm[i] = row to go in ith place
    iperm = flint_malloc(mat->r * sizeof(slong));
    for (i = 0; i < mat->r; ++i)
    {
        iperm[perm[i]] = i;
    }

    // Will do at most rows - 1 swaps
    for (i = 0; i < mat->r - 1; ++i)
    {
        // Get element to permute with current location
        for (j = iperm[i]; j < i; j = iperm[j]);
        if (i != j)
        {
            gr_sparse_vec_swap(&mat->rows[i], &mat->rows[j], ctx);
        }
    }

    flint_free(iperm);
    return GR_SUCCESS;
}
