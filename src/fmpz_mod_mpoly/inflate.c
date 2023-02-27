/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_inflate(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz * shift,
    const fmpz * stride,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int have_zero_stride;
    slong j;
    slong Abits;
    slong nvars = ctx->minfo->nvars;
    fmpz * exps;
    TMP_INIT;

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    exps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (j = 0; j < nvars; j++)
        fmpz_init(exps + j);

    /* quick and safe bound on bits required */
    mpoly_degrees_ffmpz(exps, B->exps, B->length, B->bits, ctx->minfo);
    have_zero_stride = 0;
    for (j = 0; j < nvars; j++)
    {
        have_zero_stride |= fmpz_is_zero(stride + j);
        fmpz_mul(exps + j, exps + j, stride + j);
        fmpz_add(exps + j, exps + j, shift + j);
    }
    Abits = mpoly_exp_bits_required_ffmpz(exps, ctx->minfo);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    for (j = 0; j < nvars; j++)
        fmpz_clear(exps + j);

    if (A == B)
    {
        slong NA = mpoly_words_per_exp(Abits, ctx->minfo);
        slong exps_alloc = NA*B->length;
        ulong * texps = flint_malloc(exps_alloc*sizeof(ulong));
        mpoly_monomials_inflate(texps, Abits, B->exps, B->bits, B->length,
                                                    shift, stride, ctx->minfo);
        flint_free(A->exps);
        A->exps = texps;
        A->bits = Abits;
        A->exps_alloc = exps_alloc;
    }
    else
    {
        fmpz_mod_mpoly_fit_length_reset_bits(A, B->length, Abits, ctx);
        _fmpz_vec_set(A->coeffs, B->coeffs, B->length);
        mpoly_monomials_inflate(A->exps, Abits, B->exps, B->bits, B->length,
                                                    shift, stride, ctx->minfo);
        _fmpz_mod_mpoly_set_length(A, B->length, ctx);
    }

    TMP_END;

    if (have_zero_stride || ctx->minfo->ord != ORD_LEX)
    {
        fmpz_mod_mpoly_sort_terms(A, ctx);
        if (have_zero_stride)
            fmpz_mod_mpoly_combine_like_terms(A, ctx);
    }
    return;
}
