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


int fmpz_mpoly_gcd_cofactors(
    fmpz_mpoly_t G,
    fmpz_mpoly_t Abar,
    fmpz_mpoly_t Bbar,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            fmpz_mpoly_zero(Abar, ctx);
            fmpz_mpoly_zero(Bbar, ctx);
            return 1;
        }
        fmpz_mpoly_set(G, B, ctx);
        fmpz_mpoly_zero(Abar, ctx);
        fmpz_mpoly_one(Bbar, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, G, ctx);
            fmpz_mpoly_neg(Bbar, Bbar, ctx);
        }
        return 1;
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        fmpz_mpoly_set(G, A, ctx);
        fmpz_mpoly_zero(Bbar, ctx);
        fmpz_mpoly_one(Abar, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, G, ctx);
            fmpz_mpoly_neg(Abar, Abar, ctx);
        }
        return 1;
    }

    return _fmpz_mpoly_gcd_algo(G, Abar, Bbar, A, B, ctx, MPOLY_GCD_USE_ALL);
}

