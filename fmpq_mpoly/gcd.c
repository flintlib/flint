/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


int fmpq_mpoly_gcd(fmpq_mpoly_t G, const fmpq_mpoly_t A,
                          const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
{
    int success;

    if (fmpq_mpoly_is_zero(A, ctx) && fmpq_mpoly_is_zero(B, ctx))
    {
        fmpq_mpoly_zero(G, ctx);
        return 1;
    }

    success = fmpz_mpoly_gcd(G->zpoly, A->zpoly, B->zpoly, ctx->zctx);
    if (success)
    {
        _fmpq_mpoly_make_monic_inplace(G, ctx);
    }

    return success;
}
