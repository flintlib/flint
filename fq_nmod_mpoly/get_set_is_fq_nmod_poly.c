/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


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
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong Blen = B->length;
    const mp_limb_t * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong i, N = mpoly_words_per_exp(Bbits, ctx->minfo);
    ulong k;
    int success;
    fq_nmod_t t;

    fq_nmod_init(t, ctx->fqctx);

    fq_nmod_poly_zero(A, ctx->fqctx);

    if (B->length < 1)
    {
        success = 1;
        goto cleanup;
    }

    if (Bbits <= FLINT_BITS)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - Bbits);
        slong off, shift;

        mpoly_gen_offset_shift_sp(&off, &shift, var, Bbits, ctx->minfo);

        for (i = 0; i < Blen; i++)
        {
            k = (Bexps[N*i + off] >> shift) & mask;
            n_fq_get_fq_nmod(t, Bcoeffs + d*i, ctx->fqctx);
            fq_nmod_poly_set_coeff(A, k, t, ctx->fqctx);
        }

        success = 1;
        goto cleanup;
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
            {
                success = 0;
                goto cleanup;
            }

            n_fq_get_fq_nmod(t, Bcoeffs + d*i, ctx->fqctx);
            fq_nmod_poly_set_coeff(A, k, t, ctx->fqctx);
        }

        success = 1;
        goto cleanup;
    }

cleanup:

    fq_nmod_clear(t, ctx->fqctx);
    return success;
}


void _fq_nmod_mpoly_set_fq_nmod_poly(
    fq_nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fq_nmod_struct * Bcoeffs,
    slong Blen,
    slong var,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
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

    fq_nmod_mpoly_fit_length_reset_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        if (fq_nmod_is_zero(Bcoeffs + i, ctx->fqctx))
            continue;

        FLINT_ASSERT(d*Alen < A->coeffs_alloc);
        FLINT_ASSERT(N*Alen < A->exps_alloc);
        n_fq_set_fq_nmod(A->coeffs + d*Alen, Bcoeffs + i, ctx->fqctx);
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
