/*
    Copyright (C) 2015 Anubhav Srivastava

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
fmpz_mat_concat_horizontal(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_mat_t mat2)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    GR_MUST_SUCCEED(gr_mat_concat_horizontal((gr_mat_struct *) res,
        (const gr_mat_struct *) mat1, (const gr_mat_struct *) mat2, ctx));
}

void
fmpz_mat_concat_vertical(fmpz_mat_t res, const fmpz_mat_t mat1, const fmpz_mat_t mat2)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    GR_MUST_SUCCEED(gr_mat_concat_vertical((gr_mat_struct *) res,
        (const gr_mat_struct *) mat1, (const gr_mat_struct *) mat2, ctx));
}
