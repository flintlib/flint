/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

/* exponents of B are not multiprecision */
static void _fq_zech_mpoly_evaluate_one_fq_zech_sp(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    slong var,
    const fq_zech_t val,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i, N, off, shift;
    ulong * cmpmask, * one;
    slong Blen = B->length;
    const fq_zech_struct * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t bits = B->bits;
    slong Alen;
    fq_zech_struct * Acoeffs;
    ulong * Aexps;
    ulong mask, k;
    int need_sort = 0, cmp;
    fq_zech_t pp;
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);

    TMP_START;

    fq_zech_init(pp, ctx->fqctx);

    fq_zech_mpoly_fit_length_reset_bits(A, Blen, bits, ctx);
    Acoeffs = A->coeffs;
    Aexps = A->exps;

    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_offset_shift_sp(one, &off, &shift, var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        k = (Bexps[N*i + off] >> shift) & mask;
        fq_zech_pow_ui(pp, val, k, ctx->fqctx);
        fq_zech_mul(Acoeffs + Alen, Bcoeffs + i, pp, ctx->fqctx);
        if (fq_zech_is_zero(Acoeffs + Alen, ctx->fqctx))
            continue;

        mpoly_monomial_msub(Aexps + N*Alen, Bexps + N*i, k, one, N);
        if (Alen < 1)
        {
            Alen = 1;
            continue;
        }
        cmp = mpoly_monomial_cmp(Aexps + N*(Alen - 1), Aexps + N*Alen, N, cmpmask);
        if (cmp != 0)
        {
            need_sort |= (cmp < 0);
            Alen++;
            continue;
        }
        fq_zech_add(Acoeffs + Alen - 1, Acoeffs + Alen - 1, Acoeffs + Alen, ctx->fqctx);
        Alen -= fq_zech_is_zero(Acoeffs + Alen - 1, ctx->fqctx);
    }
    A->length = Alen;

    fq_zech_clear(pp, ctx->fqctx);

    TMP_END;

    if (need_sort)
    {
        fq_zech_mpoly_sort_terms(A, ctx);
        fq_zech_mpoly_combine_like_terms(A, ctx);
    }

    FLINT_ASSERT(fq_zech_mpoly_is_canonical(A, ctx));
}


/* exponents of B are multiprecision */
static void _fq_zech_mpoly_evaluate_one_fq_zech_mp(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    slong var,
    const fq_zech_t val,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i, N, off;
    ulong * cmpmask, * one, * tmp;
    slong Blen = B->length;
    const fq_zech_struct * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t bits = B->bits;
    slong Alen;
    fq_zech_struct * Acoeffs;
    ulong * Aexps;
    fmpz_t k;
    int need_sort = 0, cmp;
    fq_zech_t pp;
    TMP_INIT;

    FLINT_ASSERT(B->bits > FLINT_BITS);
    FLINT_ASSERT((B->bits % FLINT_BITS) == 0);

    TMP_START;

    fmpz_init(k);
    fq_zech_init(pp, ctx->fqctx);

    fq_zech_mpoly_fit_length_reset_bits(A, Blen, bits, ctx);
    Acoeffs = A->coeffs;
    Aexps = A->exps;

    N = mpoly_words_per_exp_mp(bits, ctx->minfo);
    one = (ulong *) TMP_ALLOC(3*N*sizeof(ulong));
    cmpmask = one + N;
    tmp = cmpmask + N;
    off = mpoly_gen_monomial_offset_mp(one, var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        fmpz_set_ui_array(k, Bexps + N*i + off, bits/FLINT_BITS);
        fq_zech_pow(pp, val, k, ctx->fqctx);
        fq_zech_mul(Acoeffs + Alen, Bcoeffs + i, pp, ctx->fqctx);
        if (fq_zech_is_zero(Acoeffs + Alen, ctx->fqctx))
            continue;

        mpoly_monomial_mul_fmpz(tmp, one, N, k);
        mpoly_monomial_sub_mp(Aexps + N*Alen, Bexps + N*i, tmp, N);
        if (Alen < 1)
        {
            Alen = 1;
            continue;
        }
        cmp = mpoly_monomial_cmp(Aexps + N*(Alen - 1), Aexps + N*Alen, N, cmpmask);
        if (cmp != 0)
        {
            need_sort |= (cmp < 0);
            Alen++;
            continue;
        }
        fq_zech_add(Acoeffs + Alen - 1, Acoeffs + Alen - 1, Acoeffs + Alen, ctx->fqctx);
        Alen -= fq_zech_is_zero(Acoeffs + Alen - 1, ctx->fqctx);
    }
    A->length = Alen;

    fq_zech_clear(pp, ctx->fqctx);
    fmpz_clear(k);

    TMP_END;

    if (need_sort)
    {
        fq_zech_mpoly_sort_terms(A, ctx);
        fq_zech_mpoly_combine_like_terms(A, ctx);
    }

    FLINT_ASSERT(fq_zech_mpoly_is_canonical(A, ctx));
}


void fq_zech_mpoly_evaluate_one_fq_zech(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    slong var,
    const fq_zech_t val,
    const fq_zech_mpoly_ctx_t ctx)
{
    if (fq_zech_mpoly_is_zero(B, ctx))
    {
        fq_zech_mpoly_zero(A, ctx);
        return;
    }

    if (B->bits <= FLINT_BITS)
        _fq_zech_mpoly_evaluate_one_fq_zech_sp(A, B, var, val, ctx);
    else
        _fq_zech_mpoly_evaluate_one_fq_zech_mp(A, B, var, val, ctx);
}

