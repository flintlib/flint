/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/* todo: algorithm selection */
int
gr_mat_det(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    /* return gr_mat_det_lu(res, A, ctx); */
    /* return gr_mat_det_bareiss(res, A, ctx); */

    if (A->r <= 4)
        return gr_mat_det_cofactor(res, A, ctx);
    else
        return gr_mat_det_berkowitz(res, A, ctx);
}
