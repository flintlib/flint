/*
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"
#include "perm.h"

int
gr_mat_permute_cols(gr_mat_t mat, slong * perm_store, const slong * perm_act, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i;
    slong sz = ctx->sizeof_elem;
    slong * _perm_act = _perm_init(mat->c);

    /* perm_store[i] <- perm_store[perm_act[i]] */
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->c);

    /* permute each row */
    for (i = 0; i < mat->r; i++)
    {
        _perm_set(_perm_act, perm_act, mat->c);
        _gr_vec_permute_inv(GR_MAT_ENTRY(mat, i, 0, sz), _perm_act, mat->c, ctx);
    }

    _perm_clear(_perm_act);
    
    return status;
}
