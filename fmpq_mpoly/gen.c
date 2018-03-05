/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_gen(fmpq_mpoly_t poly, slong i, const fmpq_mpoly_ctx_t ctx)
{
    slong j;
    ulong * mon;
    TMP_INIT;

    TMP_START;
    fmpz_mpoly_fit_length(poly->zpoly, 1, ctx->zctx);
    fmpz_set_ui(poly->zpoly->coeffs + 0, 1);
    fmpq_one(poly->content);

    mon = (ulong *) TMP_ALLOC((ctx->zctx->minfo->nvars)*sizeof(ulong));
    for (j = 0; j < ctx->zctx->minfo->nvars; j++)
       mon[j] = (j == i);
    mpoly_set_monomial_ui(poly->zpoly->exps, mon, poly->zpoly->bits, ctx->zctx->minfo);

    TMP_END;
    _fmpz_mpoly_set_length(poly->zpoly, 1, ctx->zctx);
}
