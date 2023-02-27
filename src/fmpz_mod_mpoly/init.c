/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_init3(
    fmpz_mod_mpoly_t A,
    slong alloc,
    flint_bitcnt_t bits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    if (alloc > 0)
    {
        A->coeffs_alloc = alloc;
        A->coeffs = (fmpz *) flint_calloc(A->coeffs_alloc, sizeof(fmpz));
        A->exps_alloc = N*alloc;
        A->exps = FLINT_ARRAY_ALLOC(A->exps_alloc, ulong);
    }
    else
    {
        A->coeffs = NULL;
        A->exps = NULL;
        A->coeffs_alloc = 0;
        A->exps_alloc = 0;
    }
    A->length = 0;
    A->bits = bits;
}

void fmpz_mod_mpoly_init2(
    fmpz_mod_mpoly_t A,
    slong alloc,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(MPOLY_MIN_BITS, ctx->minfo);
    fmpz_mod_mpoly_init3(A, alloc, bits, ctx);
}
