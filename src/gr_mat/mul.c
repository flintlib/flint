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
gr_mat_mul_generic(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    return gr_mat_mul_classical(C, A, B, ctx);
}

int
gr_mat_mul(gr_mat_t res, const gr_mat_t x, const gr_mat_t y, gr_ctx_t ctx)
{
    return GR_MAT_BINARY_OP(ctx, MAT_MUL)(res, x, y, ctx);
}
