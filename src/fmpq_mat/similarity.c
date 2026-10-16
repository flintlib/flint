/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"
#include "gr.h"
#include "gr_mat.h"

void fmpq_mat_similarity(fmpq_mat_t A, slong r, fmpq_t d)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpq(ctx);
    GR_MUST_SUCCEED(gr_mat_apply_row_similarity((gr_mat_struct *) A, r, d, ctx));
}
