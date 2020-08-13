/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_mul(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                              const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx)
{
    if (fmpq_mpoly_is_zero(B, ctx) || fmpq_mpoly_is_zero(C, ctx))
    {
        fmpq_mpoly_zero(A, ctx);
        return;
    }

    fmpq_mul(A->content, B->content, C->content);
    fmpz_mpoly_mul(A->zpoly, B->zpoly, C->zpoly, ctx->zctx);
}
