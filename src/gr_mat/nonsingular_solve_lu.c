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
gr_mat_nonsingular_solve_lu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong rank, n, m, *perm;
    gr_mat_t LU;

    n = gr_mat_nrows(A, ctx);
    m = gr_mat_ncols(X, ctx);

    if (n == 0)
        return GR_SUCCESS;

    perm = _perm_init(n);
    gr_mat_init(LU, n, n, ctx);

    status = gr_mat_lu(&rank, perm, LU, A, 1, ctx);

    if (status == GR_SUCCESS && rank == n)
    {
        if (m != 0)
            status |= gr_mat_nonsingular_solve_lu_precomp(X, perm, LU, B, ctx);
    }
    else
    {
        status |= GR_DOMAIN;
    }

    gr_mat_clear(LU, ctx);
    _perm_clear(perm);

    return status;
}
