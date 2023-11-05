/*
    Copyright (C) 2010, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
arith_stirling_matrix_1u(fmpz_mat_t mat)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    GR_MUST_SUCCEED(gr_mat_stirling((gr_mat_struct *) mat, 0, ctx));
}

void
arith_stirling_matrix_1(fmpz_mat_t mat)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    GR_MUST_SUCCEED(gr_mat_stirling((gr_mat_struct *) mat, 1, ctx));
}

void
arith_stirling_matrix_2(fmpz_mat_t mat)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    GR_MUST_SUCCEED(gr_mat_stirling((gr_mat_struct *) mat, 2, ctx));
}
