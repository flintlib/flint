/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int
gr_mat_det_berkowitz(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    gr_ptr t;
    slong n;
    int status = GR_SUCCESS;

    n = A->r;

    GR_TMP_INIT_VEC(t, n + 1, ctx);

    status |= _gr_mat_charpoly_berkowitz(t, A, ctx);
    gr_swap(res, t, ctx);
    if (n % 2)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR_VEC(t, n + 1, ctx);

    return status;
}
