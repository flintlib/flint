/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_equal_fmpz(const fmpq_mpoly_t A, const fmpz_t c,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    if (fmpq_mpoly_is_zero(A, ctx))
    {
        return fmpz_is_zero(c);
    }
    else
    {
        return   fmpz_is_one(fmpq_denref(A->content))
              && fmpz_equal(fmpq_numref(A->content), c)
              && fmpz_mpoly_is_one(A->zpoly, ctx->zctx);
    }
}

int fmpq_mpoly_equal_ui(const fmpq_mpoly_t A, ulong c,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    if (fmpq_mpoly_is_zero(A, ctx))
    {
        return c == UWORD(0);
    }
    else
    {
        return   fmpz_is_one(fmpq_denref(A->content))
              && fmpz_equal_ui(fmpq_numref(A->content), c)
              && fmpz_mpoly_is_one(A->zpoly, ctx->zctx);
    }
}

int fmpq_mpoly_equal_si(const fmpq_mpoly_t A, slong c,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    if (fmpq_mpoly_is_zero(A, ctx))
    {
        return c == WORD(0);
    }
    else
    {
        return   fmpz_is_one(fmpq_denref(A->content))
              && fmpz_equal_si(fmpq_numref(A->content), c)
              && fmpz_mpoly_is_one(A->zpoly, ctx->zctx);
    }
}
