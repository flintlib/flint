/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_scalar_mul_fmpq(fmpq_mpoly_t A,
              const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_mul(A->content, B->content, c);
    if (fmpq_is_zero(A->content))
    {
        fmpz_mpoly_zero(A->zpoly, ctx->zctx);
    }
    else
    {
        fmpz_mpoly_set(A->zpoly, B->zpoly, ctx->zctx);
    }
}
