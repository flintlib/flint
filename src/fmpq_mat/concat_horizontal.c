/*
    Copyright (C) 2015 Elena Sergeicheva

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
fmpq_mat_concat_horizontal(fmpq_mat_t res, const fmpq_mat_t mat1, const fmpq_mat_t mat2)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpq(ctx);
    GR_MUST_SUCCEED(gr_mat_concat_horizontal((gr_mat_struct *) res,
        (const gr_mat_struct *) mat1, (const gr_mat_struct *) mat2, ctx));
}
