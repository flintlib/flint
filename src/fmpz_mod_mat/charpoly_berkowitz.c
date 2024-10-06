/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_mat.h"

static void _fmpz_mod_mat_charpoly_berkowitz(fmpz* cp, const fmpz_mod_mat_t mat,
                                                      const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_mat_charpoly_berkowitz(cp, (const gr_mat_struct *) mat, gr_ctx));
}

void fmpz_mod_mat_charpoly_berkowitz(fmpz_mod_poly_t cp,
                            const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    if (!fmpz_mod_mat_is_square(mat, ctx))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_mat_charpoly_berkowitz). Non-square matrix.\n");
    }

    fmpz_mod_poly_fit_length(cp, fmpz_mod_mat_nrows(mat, ctx) + 1, ctx);
    _fmpz_mod_mat_charpoly_berkowitz(cp->coeffs, mat, ctx);
    _fmpz_mod_poly_set_length(cp, fmpz_mod_mat_nrows(mat, ctx) + 1);
    _fmpz_mod_poly_normalise(cp);
}
