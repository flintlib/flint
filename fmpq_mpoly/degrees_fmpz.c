/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_degrees_fmpz(fmpz ** degs, const fmpq_mpoly_t poly,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * new_degs;
    TMP_INIT;

    TMP_START;
    new_degs = (fmpz *) TMP_ALLOC(ctx->zctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->zctx->minfo->nvars; i++)
        fmpz_init(new_degs + i);

    mpoly_degrees_fmpz(new_degs, poly->zpoly->exps, poly->zpoly->length,
                                          poly->zpoly->bits, ctx->zctx->minfo);

    for (i = 0; i < ctx->zctx->minfo->nvars; i++)
    {
        fmpz_swap(new_degs + i, degs[i]);
        fmpz_clear(new_degs + i);
    }

    TMP_END;
}
