/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_set_fmpz(fmpq_mpoly_t poly,
                                    const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
{
    if (fmpz_is_zero(c))
    {
        fmpq_mpoly_zero(poly, ctx);
        return;
    }

    fmpz_set(fmpq_numref(poly->content), c);
    fmpz_one(fmpq_denref(poly->content));
    fmpz_mpoly_one(poly->zpoly, ctx->zctx);
}

void fmpq_mpoly_set_ui(fmpq_mpoly_t poly, ulong c, const fmpq_mpoly_ctx_t ctx)
{
    if (c == UWORD(0))
    {
        fmpq_mpoly_zero(poly, ctx);
        return;
    }

    fmpz_set_ui(fmpq_numref(poly->content), c);
    fmpz_one(fmpq_denref(poly->content));
    fmpz_mpoly_one(poly->zpoly, ctx->zctx);
}

void fmpq_mpoly_set_si(fmpq_mpoly_t poly, slong c, const fmpq_mpoly_ctx_t ctx)
{
    if (c == WORD(0))
    {
        fmpq_mpoly_zero(poly, ctx);
        return;
    }

    fmpz_set_si(fmpq_numref(poly->content), c);
    fmpz_one(fmpq_denref(poly->content));
    fmpz_mpoly_one(poly->zpoly, ctx->zctx);
}

