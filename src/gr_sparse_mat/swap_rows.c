/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

int gr_lil_mat_swap_rows(gr_lil_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)
{
    /* todo: bounds checking */
    if (r < 0 || s < 0 || r >= mat->r || s >= mat->r)
    {
        return GR_DOMAIN;
    }
    if (r == s)
    {
        return GR_DOMAIN;
    }

    if (perm != NULL)
        FLINT_SWAP(slong, perm[r], perm[s]);

    if (mat->rows[r]->nnz != 0 || mat->rows[s]->nnz != 0)
    {
        gr_sparse_vec_swap(mat->rows[r], mat->rows[s], ctx);
    }

    return GR_SUCCESS;
}
