/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


void fmpz_mpoly_term_content(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t Abits;
    fmpz * minAfields, * min_degs;
    fmpz_t g;
    TMP_INIT;

    if (A->length == 0)
    {
        fmpz_mpoly_zero(M, ctx);
        return;
    }

    TMP_START;

    Abits = A->bits;

    /* get the field-wise minimum */
    minAfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_init(minAfields + i);
    mpoly_min_fields_fmpz(minAfields, A->exps, A->length, Abits, ctx->minfo);

    /* unpack to get the min exponents of each variable */
    min_degs = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(min_degs + i);
    mpoly_get_monomial_ffmpz_unpacked_ffmpz(min_degs, minAfields, ctx->minfo);

    fmpz_mpoly_fit_length_reset_bits(M, 1, Abits, ctx);
    mpoly_set_monomial_ffmpz(M->exps, min_degs, Abits, ctx->minfo);

    fmpz_init(g);
    _fmpz_vec_content(g, A->coeffs, A->length);
    fmpz_swap(M->coeffs + 0, g);
    fmpz_clear(g);

    for (i = 0; i < ctx->minfo->nfields; i++)
        fmpz_clear(minAfields + i);
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(min_degs + i);

    _fmpz_mpoly_set_length(M, 1, ctx);

    TMP_END;
}
