/*
    Copyright (C) 2019 Tommy Hofmann
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"
#include "gr.h"
#include "gr_mat.h"

slong fmpz_mod_mat_rref(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    slong rank;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_mat_rref_lu(&rank, (gr_mat_struct *) B, (const gr_mat_struct *) A, gr_ctx));
    return rank;
}
