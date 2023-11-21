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
gr_mat_det_lu(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    slong * P;
    gr_mat_t T;
    slong n, rank;
    int status;

    n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx))
        return GR_DOMAIN;

    if (n == 0)
        return gr_one(res, ctx);

    P = _perm_init(n);
    gr_mat_init(T, n, n, ctx);
    status = gr_mat_lu(&rank, P, T, A, 1, ctx);

    if (status == GR_SUCCESS)
    {
        if (rank == 0)
        {
            status |= gr_zero(res, ctx);
        }
        else
        {
            slong i, sz = ctx->sizeof_elem;

            status |= gr_set(res, GR_MAT_ENTRY(T, 0, 0, sz), ctx);
            for (i = 1; i < n; i++)
                status |= gr_mul(res, res, GR_MAT_ENTRY(T, i, i, sz), ctx);

            if (_perm_parity(P, n))
                status |= gr_neg(res, res, ctx);
        }
    }
    else
    {
        status |= GR_UNABLE;
    }

    gr_mat_clear(T, ctx);
    _perm_clear(P);

    return status;
}
