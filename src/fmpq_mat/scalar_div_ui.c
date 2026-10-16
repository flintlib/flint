/*
    Copyright (C) 2026 Lars Göttgens

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
fmpq_mat_scalar_div_ui(fmpq_mat_t res, const fmpq_mat_t mat, const ulong x)
{
    gr_ctx_t gr_ctx;
    gr_ctx_init_fmpq(gr_ctx);
    GR_MUST_SUCCEED(gr_mat_div_ui((gr_mat_struct *) res, (const gr_mat_struct *) mat, x, gr_ctx));
}
