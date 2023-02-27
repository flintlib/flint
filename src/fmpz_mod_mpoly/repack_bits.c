/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


int fmpz_mod_mpoly_repack_bits(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    flint_bitcnt_t Abits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;

    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    if (B->bits == Abits || B->length == 0)
    {
        fmpz_mod_mpoly_set(A, B, ctx);
        return 1;
    }

    if (A == B)
        return fmpz_mod_mpoly_repack_bits_inplace(A, Abits, ctx);

    fmpz_mod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);

    success = mpoly_repack_monomials(A->exps, Abits, B->exps, B->bits,
                                                        B->length, ctx->minfo);
    if (success)
    {
        _fmpz_vec_set(A->coeffs, B->coeffs, B->length);
        A->length = B->length;
    }
    else
    {
        A->length = 0;
    }

    return success;
}


int fmpz_mod_mpoly_repack_bits_inplace(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    ulong * texps;
    slong talloc;
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);

    if (A->bits == Abits)
    {
        return 1;
    }

    if (A->length < 1)
    {
        A->bits = Abits;
        return 1;
    }

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    talloc = N*A->length;
    texps = FLINT_ARRAY_ALLOC(talloc, ulong);
    success = mpoly_repack_monomials(texps, Abits,
                                      A->exps, A->bits, A->length, ctx->minfo);
    A->bits = Abits;
    if (success)
    {
        flint_free(A->exps);
        A->exps = texps;
        A->exps_alloc = talloc;
    }
    else
    {
        flint_free(texps);
        A->length = 0;
    }
    return success;
}
