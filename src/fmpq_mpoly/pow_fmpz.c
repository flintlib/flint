/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_pow_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                    const fmpz_t k, const fmpq_mpoly_ctx_t ctx)
{
    int success;

    success = fmpq_pow_fmpz(A->content, B->content, k);
    success = success && fmpz_mpoly_pow_fmpz(A->zpoly, B->zpoly, k, ctx->zctx);

    return success;
}
