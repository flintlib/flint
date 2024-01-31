/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "gr.h"
#include "gr_mat.h"

void fmpz_mod_mat_det(fmpz_t res, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    slong n = mat->r;

    if (!fmpz_mod_mat_is_square(mat, ctx))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_mat_charpoly_berkowitz). Non-square matrix.\n");
    }

    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);

    if (n <= 4)
    {
        GR_MUST_SUCCEED(gr_mat_det_cofactor(res, (const gr_mat_struct *) mat, gr_ctx));
    }
    else
    {
        if (gr_mat_det_lu(res, (const gr_mat_struct *) mat, gr_ctx) != GR_SUCCESS)
        {
            /* Fall back on division-free algorithm if we encountered an impossible inverse */
            /* Could try something else here: faddeev_bsgs (O(n^3.5)) or Howell form. */
            GR_MUST_SUCCEED(gr_mat_det_berkowitz(res, (const gr_mat_struct *) mat, gr_ctx));
        }
    }
}