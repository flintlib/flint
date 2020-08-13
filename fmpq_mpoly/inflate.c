/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_inflate(fmpq_mpoly_t A, const fmpq_mpoly_t B,
           const fmpz * shift, const fmpz * stride, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_mpoly_inflate(A->zpoly, A->zpoly, shift, stride, ctx->zctx);
    fmpq_set(A->content, B->content);
    fmpq_mpoly_reduce(A, ctx);
}
