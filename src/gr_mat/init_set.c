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
gr_mat_init_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong r, c;

    r = gr_mat_nrows(mat, ctx);
    c = gr_mat_ncols(mat, ctx);

    gr_mat_init(res, r, c, ctx);
    return gr_mat_set(res, mat, ctx);

}
