/*
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
fmpz_mod_mat_pow_ui(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, ulong exp,
                    const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;

    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_mat_pow_ui
                    ((gr_mat_struct *) B, (const gr_mat_struct *) A, exp,
                     gr_ctx));
}
