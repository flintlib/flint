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
gr_mat_rref_den(slong * res_rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx)
{
    return gr_mat_rref_den_fflu(res_rank, R, den, A, ctx);
}
