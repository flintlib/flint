/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


int fmpq_mpoly_gcd_threaded(fmpq_mpoly_t G, const fmpq_mpoly_t A,
         const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx, slong thread_limit)
{
    int success;

    success = fmpz_mpoly_gcd_threaded(G->zpoly, A->zpoly, B->zpoly,
                                                      ctx->zctx, thread_limit);
    if (!success)
        return 0;

    if (G->zpoly->length > 0)
    {
        fmpz_one(fmpq_numref(G->content));
        fmpz_set(fmpq_denref(G->content), G->zpoly->coeffs + 0);
    }
    else
    {
        fmpq_zero(G->content);
    }

    return 1;
}

int fmpq_mpoly_gcd(fmpq_mpoly_t G, const fmpq_mpoly_t A,
                              const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
{
    return fmpq_mpoly_gcd_threaded(G, A, B, ctx, MPOLY_DEFAULT_THREAD_LIMIT);
}
