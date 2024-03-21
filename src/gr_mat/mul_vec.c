
/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"
#include "gr_vec.h"

int gr_mat_mul_vec(gr_ptr v, const gr_mat_t A, gr_srcptr u, gr_ctx_t ctx)
{
    slong r, c, row, sz;
    gr_ptr w;
    int status;

    sz = ctx->sizeof_elem;
    r = gr_mat_nrows(A, ctx);
    c = gr_mat_ncols(A, ctx);

    if (u == v)
    {
        GR_TMP_INIT_VEC(w, r, ctx);
        _gr_vec_init(w, r, ctx);
        status = gr_mat_mul_vec(w, A, u, ctx);
        _gr_vec_swap(v, w, r, ctx);
        GR_TMP_CLEAR_VEC(w, r, ctx);
        return status;
    }
    status = _gr_vec_zero(v, r, ctx);
    for (row = 0; row < r; ++row)
        status |= _gr_vec_dot(GR_ENTRY(v, row, sz), GR_ENTRY(v, row, sz), 0, A->rows[row], u, c, ctx);
    return status;
}
