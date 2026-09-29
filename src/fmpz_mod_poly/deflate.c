/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_poly.h"

void fmpz_mod_poly_deflate(fmpz_mod_poly_t result, const fmpz_mod_poly_t input,
                                     ulong deflation, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;

    if (deflation == 0)
        flint_throw(FLINT_DIVZERO, "fmpz_mod_poly_deflate");

    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(gr_poly_deflate((gr_poly_struct *) result,
                        (const gr_poly_struct *) input, deflation, gr_ctx));
}
