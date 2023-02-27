/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_make_monic(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t c;

    if (B->length < 1)
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_make_monic: polynomial is zero");

    fmpz_init(c);
    fmpz_mod_inv(c, B->coeffs + 0, ctx->ffinfo);
    fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(A, B, c, ctx);
    fmpz_clear(c);
}
