/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_set_fmpz(fmpq_mpoly_t A,
                                    const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_set(fmpq_numref(A->content), c);
    fmpz_one(fmpq_denref(A->content));
    if (fmpz_is_zero(c))
    {
        fmpz_mpoly_zero(A->zpoly, ctx->zctx);
    }
    else
    {
        fmpz_mpoly_one(A->zpoly, ctx->zctx);
    }
}

void fmpq_mpoly_set_ui(fmpq_mpoly_t A, ulong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_set_ui(fmpq_numref(A->content), c);
    fmpz_one(fmpq_denref(A->content));
    if (c == UWORD(0))
    {
        fmpz_mpoly_zero(A->zpoly, ctx->zctx);
    }
    else
    {
        fmpz_mpoly_one(A->zpoly, ctx->zctx);
    }
}

void fmpq_mpoly_set_si(fmpq_mpoly_t A, slong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_set_si(fmpq_numref(A->content), c);
    fmpz_one(fmpq_denref(A->content));
    if (c == WORD(0))
    {
        fmpz_mpoly_zero(A->zpoly, ctx->zctx);
    }
    else
    {
        fmpz_mpoly_one(A->zpoly, ctx->zctx);
    }
}

