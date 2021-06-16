/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int nmod_mpoly_is_nmod_poly(
    const nmod_mpoly_t A,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    return mpoly_is_poly(A->exps, A->length, A->bits, var, ctx->minfo);
}

int nmod_mpoly_get_nmod_poly(
    nmod_poly_t A,
    const nmod_mpoly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    A->mod = ctx->mod;
    return nmod_mpoly_get_n_poly(evil_cast_nmod_poly_to_n_poly(A), B, var, ctx);
}

int nmod_mpoly_get_n_poly(
    n_poly_t A,
    const nmod_mpoly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    slong Blen = B->length;
    const mp_limb_t * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t Bbits = B->bits;
    slong i, N = mpoly_words_per_exp(Bbits, ctx->minfo);
    ulong k;

    n_poly_zero(A);

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
            n_poly_set_coeff(A, k, Bcoeffs[i]);
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

            n_poly_set_coeff(A, k, Bcoeffs[i]);
        }
        return 1;
    }
}

void _nmod_mpoly_set_nmod_poly(
    nmod_mpoly_t A,
    flint_bitcnt_t Abits,
    const mp_limb_t * Bcoeffs,
    slong Blen,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    slong i, Alen;
    ulong * genexp;
    TMP_INIT;

    TMP_START;

    genexp = TMP_ARRAY_ALLOC(N, ulong);
    if (Abits <= FLINT_BITS)
        mpoly_gen_monomial_sp(genexp, var, Abits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(genexp, var, Abits, ctx->minfo);

    Alen = 2;
    for (i = 0; i < Blen; i++)
        Alen += (Bcoeffs[i] != 0);

    nmod_mpoly_fit_length_reset_bits(A, Alen, Abits, ctx);

    Alen = 0;
    for (i = Blen - 1; i >= 0; i--)
    {
        if (Bcoeffs[i] == 0)
            continue;

        FLINT_ASSERT(Alen < A->coeffs_alloc);
        A->coeffs[Alen] = Bcoeffs[i];
        if (Abits <= FLINT_BITS)
            mpoly_monomial_mul_ui(A->exps + N*Alen, genexp, N, i);
        else
            mpoly_monomial_mul_ui_mp(A->exps + N*Alen, genexp, N, i);
        Alen++;
    }
    A->length = Alen;

    TMP_END;
}

void nmod_mpoly_set_n_poly_mod(
    nmod_mpoly_t A,
    const n_poly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;

    if (B->length < 1)
    {
        nmod_mpoly_zero(A, ctx);
        return;
    }

    bits = mpoly_gen_pow_exp_bits_required(var, B->length - 1, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);
    _nmod_mpoly_set_nmod_poly(A, bits, B->coeffs, B->length, var, ctx);
}

void nmod_mpoly_set_nmod_poly(
    nmod_mpoly_t A,
    const nmod_poly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_set_n_poly_mod(A, evil_const_cast_nmod_poly_to_n_poly(B), var, ctx);
}
