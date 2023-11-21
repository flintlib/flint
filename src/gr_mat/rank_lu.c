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
gr_mat_rank_lu(slong * rank, const gr_mat_t A, gr_ctx_t ctx)
{
    slong n, m;
    slong * P;
    int status;
    gr_mat_t T;

    n = gr_mat_nrows(A, ctx);
    m = gr_mat_ncols(A, ctx);

    if (n == 0 || m == 0)
    {
        *rank = 0;
        return GR_SUCCESS;
    }
    else
    {
        gr_mat_init(T, n, m, ctx);
        P = _perm_init(n);

        status = gr_mat_lu(rank, P, T, A, 0, ctx);

        gr_mat_clear(T, ctx);
        _perm_clear(P);

        if (status != GR_SUCCESS)
            status |= GR_UNABLE;

        return status;
    }
}
