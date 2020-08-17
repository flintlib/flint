/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"
#include "nmod_mpoly_factor.h"


int fq_nmod_mpoly_is_fq_nmod_poly(
    const fq_nmod_mpoly_t A,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    return mpoly_is_poly(A->exps, A->length, A->bits, var, ctx->minfo);
}


int fq_nmod_mpoly_get_fq_nmod_poly(
    fq_nmod_poly_t A,
    const fq_nmod_mpoly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong Blen = B->length;
    const fq_nmod_struct * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong i, N = mpoly_words_per_exp(Bbits, ctx->minfo);
    ulong k;

    fq_nmod_poly_zero(A, ctx->fqctx);

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
            fq_nmod_poly_set_coeff(A, k, Bcoeffs + i, ctx->fqctx);
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

            fq_nmod_poly_set_coeff(A, k, Bcoeffs + i, ctx->fqctx);
        }
        return 1;
    }
}

void fq_nmod_mpoly_fit_length_set_bits(
    fq_nmod_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong new_alloc;

    FLINT_ASSERT(len >= 0);

    if (A->alloc < len)
    {
        new_alloc = FLINT_MAX(len, 2*A->alloc);
        if (A->alloc > 0)
        {
            A->coeffs = (fq_nmod_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(fq_nmod_struct));
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
        }
        else
        {
            A->coeffs = (fq_nmod_struct *) flint_malloc(new_alloc*sizeof(fq_nmod_struct));
            A->exps   = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
        }
        for (i = A->alloc; i < new_alloc; i++)
            fq_nmod_init(A->coeffs + i, ctx->fqctx);
        A->alloc = new_alloc;
    }
    else if (A->bits < bits)
    {
        if (A->alloc > 0)
            A->exps = (ulong *) flint_realloc(A->exps, A->alloc*N*sizeof(ulong));
    }

    A->bits = bits;
}

void _fq_nmod_mpoly_set_fq_nmod_poly(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_struct * Bcoeffs,
    slong Blen,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
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
        Alen += !fq_nmod_is_zero(Bcoeffs + i, ctx->fqctx);

    fq_nmod_mpoly_fit_length_set_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        if (fq_nmod_is_zero(Bcoeffs + i, ctx->fqctx))
            continue;

        FLINT_ASSERT(Alen < A->alloc);
        fq_nmod_set(A->coeffs + Alen, Bcoeffs + i, ctx->fqctx);
        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    A->length = Alen;

    TMP_END;
}

void fq_nmod_mpoly_set_fq_nmod_poly(
    fq_nmod_mpoly_t A,
    const fq_nmod_poly_t B,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;

    if (B->length < 1)
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    bits = mpoly_gen_pow_exp_bits_required(var, B->length - 1, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    _fq_nmod_mpoly_set_fq_nmod_poly(A, bits, B->coeffs, B->length, var, ctx);
}
