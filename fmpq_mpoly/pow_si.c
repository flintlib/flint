/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_pow_si(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                           slong k, const fmpq_mpoly_ctx_t ctx)
{
    if (k < 0)
    {
        flint_throw(FLINT_ERROR, "Negative power in fmpq_mpoly_pow_si");
    }

    fmpq_pow_si(A->content, B->content, k);
    fmpz_mpoly_pow_fps(A->zpoly, B->zpoly, k, ctx->zctx);
    return;
}
