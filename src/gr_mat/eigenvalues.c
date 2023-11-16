/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int
gr_mat_eigenvalues(gr_vec_t lambda, gr_vec_t mult, const gr_mat_t mat, int flags, gr_ctx_t ctx)
{
    int status;
    gr_poly_t cp;
    gr_poly_init(cp, ctx);
    status = gr_mat_charpoly(cp, mat, ctx);
    if (status == GR_SUCCESS)
        status = gr_poly_roots(lambda, mult, cp, 0, ctx);
    gr_poly_clear(cp, ctx);
    return status;
}

int
gr_mat_eigenvalues_other(gr_vec_t lambda, gr_vec_t mult, const gr_mat_t mat, gr_ctx_t mat_ctx, int flags, gr_ctx_t ctx)
{
    int status;
    gr_poly_t cp;
    gr_poly_init(cp, mat_ctx);
    status = gr_mat_charpoly(cp, mat, mat_ctx);
    if (status == GR_SUCCESS)
        status = gr_poly_roots_other(lambda, mult, cp, mat_ctx, 0, ctx);
    gr_poly_clear(cp, mat_ctx);
    return status;
}
