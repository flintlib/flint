/*
    Copyright (C) 2018,2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

static void _make_monic(fmpq_mpoly_t G)
{
    if (G->zpoly->length > 0)
    {
        fmpz_one(fmpq_numref(G->content));
        fmpz_set(fmpq_denref(G->content), G->zpoly->coeffs + 0);
    }
    else
    {
        fmpq_zero(G->content);
    }
}

#define FMPQ_MPOLY_GCD_EXT(ext)                                             \
int fmpq_mpoly_##ext(fmpq_mpoly_t G, const fmpq_mpoly_t A,                  \
                         const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)  \
{                                                                           \
    int success;                                                            \
    success = fmpz_mpoly_##ext(G->zpoly, A->zpoly, B->zpoly, ctx->zctx);    \
    if (success)                                                            \
        _make_monic(G);                                                     \
    return success;                                                         \
}

FMPQ_MPOLY_GCD_EXT(gcd)
FMPQ_MPOLY_GCD_EXT(gcd_brown)
FMPQ_MPOLY_GCD_EXT(gcd_hensel)
FMPQ_MPOLY_GCD_EXT(gcd_subresultant)
FMPQ_MPOLY_GCD_EXT(gcd_zippel)
FMPQ_MPOLY_GCD_EXT(gcd_zippel2)

