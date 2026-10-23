/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_poly.h"

void fmpz_mod_poly_make_monic(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                                      const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;

    if (poly->length == 0)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);

    if (gr_poly_make_monic((gr_poly_struct *) res,
            (const gr_poly_struct *) poly, gr_ctx) != GR_SUCCESS)
        flint_throw(FLINT_IMPINV, "Exception in fmpz_mod_poly_make_monic: "
                                  "leading coefficient is not invertible.\n");
}
