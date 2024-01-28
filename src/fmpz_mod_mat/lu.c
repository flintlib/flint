/*
    Copyright (C) 2011, 2024 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"
#include "fmpz_mod_mat.h"

slong fmpz_mod_mat_lu(slong * P, fmpz_mod_mat_t A, int rank_check, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    slong rank;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_mat_lu(&rank, P, (gr_mat_struct *) A, (const gr_mat_struct *) A, rank_check, gr_ctx));
    return rank;
}
