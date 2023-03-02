/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"

int fmpz_mpoly_gcd_brown(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx) || fmpz_mpoly_is_zero(B, ctx))
        return fmpz_mpoly_gcd(G, A, B, ctx);

    return _fmpz_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_BROWN);
}

