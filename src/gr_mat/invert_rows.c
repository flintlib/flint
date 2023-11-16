/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int gr_mat_invert_rows(gr_mat_t mat, slong * perm, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong r = gr_mat_nrows(mat, ctx);
    slong i;

    for (i = 0; i < r / 2; i++)
        status |= gr_mat_swap_rows(mat, perm, i, r - i - 1, ctx);

    return status;
}
