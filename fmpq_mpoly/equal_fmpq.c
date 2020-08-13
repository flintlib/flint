/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_equal_fmpq(const fmpq_mpoly_t A, const fmpq_t c,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    if (fmpq_mpoly_is_zero(A, ctx))
    {
        return fmpq_is_zero(c);
    }
    else
    {
        return   fmpq_equal(A->content, c)
              && fmpz_mpoly_is_one(A->zpoly, ctx->zctx);
    }
}
