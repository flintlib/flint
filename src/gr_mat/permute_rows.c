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
gr_mat_permute_rows(gr_mat_t mat, slong * perm_store, const slong * perm_act, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i;
    slong sz = ctx->sizeof_elem;
    gr_ptr mat_tmp;

    GR_TMP_INIT_VEC(mat_tmp, mat->r * mat->c, ctx);

    /* perm_store[i] <- perm_store[perm_act[i]] */
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->r);

    /* rows[i] <- rows[perm_act[i]]  */
    for (i = 0; i < mat->r; i++)
        status |= _gr_vec_set(GR_ENTRY(mat_tmp, i * mat->c, sz), GR_MAT_ENTRY(mat, perm_act[i], 0, sz), mat->c, ctx);

    for (i = 0; i < mat->r; i++)
        status |= _gr_vec_set(GR_MAT_ENTRY(mat, i, 0, sz), GR_ENTRY(mat_tmp, i * mat->c, sz), mat->c, ctx);

    GR_TMP_CLEAR_VEC(mat_tmp, mat->r * mat->c, ctx);
    
    return status;
}
