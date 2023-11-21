/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_submul_scalar(gr_mat_t res, const gr_mat_t mat, gr_srcptr x, gr_ctx_t ctx)
{
    slong i, r, c;
    int status = GR_SUCCESS;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (c != 0)
        for (i = 0; i < r; i++)
            status |= _gr_vec_submul_scalar(res->rows[i], mat->rows[i], c, x, ctx);

    return status;
}
