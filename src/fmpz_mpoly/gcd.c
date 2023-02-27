/*
    Copyright (C) 2018-2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"

int fmpz_mpoly_gcd(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
            fmpz_mpoly_zero(G, ctx);
        else if (fmpz_sgn(B->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, B, ctx);
        else
            fmpz_mpoly_set(G, B, ctx);

        return 1;
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, A, ctx);
        else
            fmpz_mpoly_set(G, A, ctx);

        return 1;
    }

    return _fmpz_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ALL);
}

