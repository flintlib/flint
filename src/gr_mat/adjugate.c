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
gr_mat_adjugate(gr_mat_t adj, gr_ptr det, const gr_mat_t A, gr_ctx_t ctx)
{
    if (gr_mat_nrows(A, ctx) <= 5)
        return gr_mat_adjugate_cofactor(adj, det, A, ctx);
    else
        return gr_mat_adjugate_charpoly(adj, det, A, ctx);
}
