/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


int fmpq_mpoly_repack_bits(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                 flint_bitcnt_t Abits, const fmpq_mpoly_ctx_t ctx)
{
    int success;

    success = fmpz_mpoly_repack_bits(A->zpoly, B->zpoly, Abits, ctx->zctx);
    if (success)
    {
        fmpq_set(A->content, B->content);
    }

    return success;
}
