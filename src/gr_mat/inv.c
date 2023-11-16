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
gr_mat_inv(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    gr_mat_t T;
    int status = GR_SUCCESS;
    slong n = mat->r;

    if (n != mat->c)
        return GR_DOMAIN;

    if (n == 0)
        return GR_SUCCESS;

    if (n == 1)
        return gr_inv(res->rows[0], mat->rows[0], ctx);

    gr_mat_init(T, n, n, ctx);

    status |= gr_mat_one(T, ctx);
    status |= gr_mat_nonsingular_solve(res, mat, T, ctx);

    gr_mat_clear(T, ctx);

    return status;
}
