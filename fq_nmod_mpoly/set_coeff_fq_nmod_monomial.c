/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_set_coeff_fq_nmod_monomial(fq_nmod_mpoly_t A,
                               const fq_nmod_t c, const fq_nmod_mpoly_t M,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    fmpz * texps;
    TMP_INIT;

    if (M->length != WORD(1))
    {
        flint_throw(FLINT_ERROR,
                 "M not monomial in fq_nmod_mpoly_set_coeff_fq_nmod_monomial");
    }

    TMP_START;
    texps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (i = 0; i < nvars; i++)
        fmpz_init(texps + i);

    mpoly_get_monomial_ffmpz(texps, M->exps + 0, M->bits, ctx->minfo);
    _fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(A, c, texps, ctx);

    for (i = 0; i < nvars; i++)
        fmpz_clear(texps + i);

    TMP_END;
    return;
}
