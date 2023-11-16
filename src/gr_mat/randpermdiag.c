/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "gr_mat.h"

int
gr_mat_randpermdiag(int * _parity, gr_mat_t mat, flint_rand_t state, gr_ptr diag, slong n, gr_ctx_t ctx)
{
    int parity;
    slong i;
    slong *rows;
    slong *cols;
    int status = GR_SUCCESS;

    if (n > mat->r || n > mat->c)
        return GR_DOMAIN;

    rows = _perm_init(mat->r);
    cols = _perm_init(mat->c);

    parity = _perm_randtest(rows, mat->r, state);
    parity ^= _perm_randtest(cols, mat->c, state);

    status |= gr_mat_zero(mat, ctx);
    for (i = 0; i < n; i++)
        status |= gr_set(gr_mat_entry_ptr(mat, rows[i], cols[i], ctx), GR_ENTRY(diag, i, ctx->sizeof_elem), ctx);

    _perm_clear(rows);
    _perm_clear(cols);

    *_parity = parity;
    return status;
}
