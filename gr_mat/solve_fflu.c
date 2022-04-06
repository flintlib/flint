/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/perm.h"
#include "gr_mat.h"

int
gr_mat_solve_fflu(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    int status;
    gr_ptr den;
    slong m;

    GR_TMP_START;
    GR_TMP_INIT1(den, ctx);

    status = gr_mat_solve_den_fflu(X, den, A, B, ctx);

    m = gr_mat_ncols(X, ctx);

    if (status == GR_SUCCESS)
    {
        if (m != 0)
            status |= gr_mat_div_scalar(X, X, den, ctx);
    }

    GR_TMP_CLEAR1(den, ctx);
    GR_TMP_END;

    return status;
}
