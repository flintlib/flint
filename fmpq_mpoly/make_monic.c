/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_make_monic(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    if (B->zpoly->length == 0)
    {
        flint_throw(FLINT_ERROR, "Zero polynomial in fmpq_mpoly_make_monic");
    }

    fmpz_one(fmpq_numref(A->content));
    fmpz_set(fmpq_denref(A->content), B->zpoly->coeffs + 0);

    if (A != B)
    {
        fmpz_mpoly_set(A->zpoly, B->zpoly, ctx->zctx);
    }
}

/* this version must ignore A->content */
void _fmpq_mpoly_make_monic_inplace(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->zpoly->length > 0);
    FLINT_ASSERT(fmpz_sgn(A->zpoly->coeffs + 0) > 0);

    fmpz_one(fmpq_numref(A->content));
    fmpz_set(fmpq_denref(A->content), A->zpoly->coeffs + 0);
}
