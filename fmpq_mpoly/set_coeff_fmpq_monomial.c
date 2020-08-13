/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_set_coeff_fmpq_monomial(fmpq_mpoly_t poly1, const fmpq_t c,
                         const fmpq_mpoly_t poly2, const fmpq_mpoly_ctx_t ctx)
{


    slong i, nvars = ctx->zctx->minfo->nvars;
    fmpz * texps;
    TMP_INIT;

    if (poly2->zpoly->length != WORD(1))
    {
        flint_throw(FLINT_ERROR, "poly2 not monomial in fmpz_mpoly_set_coeff_fmpz_monomial");
    }

    TMP_START;
    texps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
        fmpz_init(texps + i);

    mpoly_get_monomial_ffmpz(texps, poly2->zpoly->exps + 0,
                                         poly2->zpoly->bits, ctx->zctx->minfo);
    _fmpq_mpoly_set_coeff_fmpq_fmpz(poly1, c, texps, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(texps + i);

    TMP_END;
    return;
}
