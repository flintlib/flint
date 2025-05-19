/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
fmpz_mat_transpose(fmpz_mat_t B, const fmpz_mat_t A)
{
    gr_ctx_t gr_ctx;
    gr_ctx_init_fmpz(gr_ctx);
    GR_MUST_SUCCEED(gr_mat_transpose((gr_mat_struct *)B, (const gr_mat_struct *) A, gr_ctx));
}
