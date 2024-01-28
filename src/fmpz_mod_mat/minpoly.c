/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_mat.h"

void fmpz_mod_mat_minpoly(fmpz_mod_poly_t p, const fmpz_mod_mat_t X,
                                                      const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    slong n = fmpz_mod_mat_nrows(X, ctx);

    if (n != fmpz_mod_mat_ncols(X, ctx))
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_mat_minpoly). Non-square matrix.\n");

    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_mat_minpoly_field((gr_poly_struct *) p, (const gr_mat_struct *) X, gr_ctx));
}
