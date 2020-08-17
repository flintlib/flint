/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "nmod_mpoly_factor.h"


int fmpz_mpoly_is_fmpz_poly(
    const fmpz_mpoly_t A,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    return mpoly_is_poly(A->exps, A->length, A->bits, var, ctx->minfo);
}


int fmpz_mpoly_get_fmpz_poly(
    fmpz_poly_t A,
    const fmpz_mpoly_t B,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    slong Blen = B->length;
    const fmpz * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong i, N = mpoly_words_per_exp(Bbits, ctx->minfo);
    ulong k;

    fmpz_poly_zero(A);

    if (B->length < 1)
        return 1;

    if (Bbits <= FLINT_BITS)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - Bbits);
        slong off, shift;

        mpoly_gen_offset_shift_sp(&off, &shift, var, Bbits, ctx->minfo);

        for (i = 0; i < Blen; i++)
        {
            k = (Bexps[N*i + off] >> shift) & mask;
            fmpz_poly_set_coeff_fmpz(A, k, Bcoeffs + i);
        }
        return 1;
    }
    else
    {
        slong j, off;
        ulong check, wpf = Bbits/FLINT_BITS;

        off = mpoly_gen_offset_mp(var, Bbits, ctx->minfo);

        for (i = 0; i < Blen; i++)
        {
            k = Bexps[N*i + off + 0];
            check = 0;
            for (j = 1; j < wpf; j++)
                check |= Bexps[N*i + off + j];

            if (check != 0 || (slong) k < 0)
                return 0;

            fmpz_poly_set_coeff_fmpz(A, k, Bcoeffs + i);
        }
        return 1;
    }
}

void fmpz_mpoly_fit_length_set_bits(
    fmpz_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong new_alloc;

    FLINT_ASSERT(len >= 0);

    if (A->alloc < len)
    {
        new_alloc = FLINT_MAX(len, 2*A->alloc);
        if (A->alloc > 0)
        {
            A->coeffs = (fmpz *) flint_realloc(A->coeffs, new_alloc*sizeof(fmpz));
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
            memset(A->coeffs + A->alloc, 0, (new_alloc - A->alloc)*sizeof(fmpz));
        }
        else
        {
            A->coeffs = (fmpz *) flint_calloc(new_alloc, sizeof(fmpz));
            A->exps   = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
        }
        A->alloc = new_alloc;
    }
    else if (A->bits < bits)
    {
        if (A->alloc > 0)
            A->exps = (ulong *) flint_realloc(A->exps, A->alloc*N*sizeof(ulong));
    }

    A->bits = bits;
}

void _fmpz_mpoly_set_fmpz_poly(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz * Bcoeffs,
    slong Blen,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    slong i, Alen;
    ulong * genexp;
    TMP_INIT;

    TMP_START;

    genexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    if (Abits <= FLINT_BITS)
        mpoly_gen_monomial_sp(genexp, var, Abits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(genexp, var, Abits, ctx->minfo);

    Alen = 2;
    for (i = 0; i < Blen; i++)
        Alen += (Bcoeffs[i] != 0);

    fmpz_mpoly_fit_length_set_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        if (fmpz_is_zero(Bcoeffs + i))
            continue;

        FLINT_ASSERT(Alen < A->alloc);
        fmpz_set(A->coeffs + Alen, Bcoeffs + i);
        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}

void fmpz_mpoly_set_fmpz_poly(
    fmpz_mpoly_t A,
    const fmpz_poly_t B,
    slong v,
    const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;

    if (B->length < 1)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    bits = mpoly_gen_pow_exp_bits_required(v, B->length - 1, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    _fmpz_mpoly_set_fmpz_poly(A, bits, B->coeffs, B->length, v, ctx);
}
