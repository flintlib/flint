/*
    Copyright (C) 2018, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

/* exponents of B are not multiprecision */
static int _fmpq_mpoly_evaluate_one_fmpq_sp(
    fmpq_mpoly_t A,
    const fmpq_mpoly_t B,
    slong var,
    fmpz_pow_cache_t num_cache,
    fmpz_pow_cache_t den_cache,
    ulong deg,
    const fmpq_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i, N, off, shift;
    ulong * cmpmask, * one;
    slong Blen = B->zpoly->length;
    const fmpz * Bcoeffs = B->zpoly->coeffs;
    const ulong * Bexps = B->zpoly->exps;
    flint_bitcnt_t bits = B->zpoly->bits;
    slong Alen;
    fmpz * Acoeffs;
    ulong * Aexps;
    ulong mask, k;
    int need_sort = 0, cmp;
    TMP_INIT;

    FLINT_ASSERT(bits <= FLINT_BITS);

    TMP_START;

    if (A != B)
        fmpz_mpoly_fit_length_reset_bits(A->zpoly, Blen, bits, ctx->zctx);

    Acoeffs = A->zpoly->coeffs;
    Aexps = A->zpoly->exps;

    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    N = mpoly_words_per_exp_sp(bits, ctx->zctx->minfo);
    cmpmask = TMP_ARRAY_ALLOC(2*N, ulong);
    one = cmpmask + N;
    mpoly_gen_monomial_offset_shift_sp(one, &off, &shift, var, bits, ctx->zctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->zctx->minfo);

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        k = (Bexps[N*i + off] >> shift) & mask;
        success = fmpz_pow_cache_mulpow_ui(Acoeffs + Alen, Bcoeffs + i,
                                                              k, num_cache) &&
                  fmpz_pow_cache_mulpow_ui(Acoeffs + Alen, Acoeffs + Alen,
                                                           deg - k, den_cache);
        if (!success)
            break;

        if (fmpz_is_zero(Acoeffs + Alen))
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
        fmpz_add(Acoeffs + Alen - 1, Acoeffs + Alen - 1, Acoeffs + Alen);
        Alen -= fmpz_is_zero(Acoeffs + Alen - 1);
    }

    /* from the fmpz_add: at most two junk coeffs past length */
    for (i = Alen; i < Alen + 2 && i < A->zpoly->alloc; i++)
        _fmpz_demote(Acoeffs + i);

    _fmpz_mpoly_set_length(A->zpoly, Alen, ctx->zctx);

    TMP_END;

    if (success)
    {
        fmpz_set(fmpq_numref(A->content), fmpq_numref(B->content));
        success = fmpz_pow_cache_mulpow_ui(fmpq_denref(A->content),
                                      fmpq_denref(B->content), deg, den_cache);
        fmpq_canonicalise(A->content);
    }
    else
    {
        fmpq_zero(A->content);
    }

    if (need_sort)
    {
        fmpz_mpoly_sort_terms(A->zpoly, ctx->zctx);
        fmpz_mpoly_combine_like_terms(A->zpoly, ctx->zctx);
    }

    fmpq_mpoly_reduce(A, ctx);

    FLINT_ASSERT(fmpq_mpoly_is_canonical(A, ctx));

    return success;
}

/* exponents of B are multiprecision */
static int _fmpq_mpoly_evaluate_one_fmpq_mp(
    fmpq_mpoly_t A,
    const fmpq_mpoly_t B,
    slong var,
    fmpz_pow_cache_t num_cache,
    fmpz_pow_cache_t den_cache,
    const fmpz_t deg,
    const fmpq_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i, N, off;
    ulong * cmpmask, * one, * tmp;
    slong Blen = B->zpoly->length;
    const fmpz * Bcoeffs = B->zpoly->coeffs;
    const ulong * Bexps = B->zpoly->exps;
    flint_bitcnt_t bits = B->zpoly->bits;
    slong Alen;
    fmpz * Acoeffs;
    ulong * Aexps;
    fmpz_t k, degmk;
    int need_sort = 0, cmp;
    TMP_INIT;

    FLINT_ASSERT((bits % FLINT_BITS) == 0);

    TMP_START;

    fmpz_init(k);
    fmpz_init(degmk);
    
    if (A != B)
        fmpz_mpoly_fit_length_reset_bits(A->zpoly, Blen, bits, ctx->zctx);

    Acoeffs = A->zpoly->coeffs;
    Aexps = A->zpoly->exps;

    N = mpoly_words_per_exp(bits, ctx->zctx->minfo);
    one = (ulong *) TMP_ALLOC(3*N*sizeof(ulong));
    cmpmask = one + N;
    tmp = cmpmask + N;
    off = mpoly_gen_monomial_offset_mp(one, var, bits, ctx->zctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->zctx->minfo);

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        fmpz_set_ui_array(k, Bexps + N*i + off, bits/FLINT_BITS);
        fmpz_sub(degmk, deg, k);
        success = fmpz_pow_cache_mulpow_fmpz(Acoeffs + Alen, Bcoeffs + i,
                                                              k, num_cache) &&
                  fmpz_pow_cache_mulpow_fmpz(Acoeffs + Alen, Acoeffs + Alen,
                                                             degmk, den_cache);
        if (!success)
            break;

        if (fmpz_is_zero(Acoeffs + Alen))
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
        fmpz_add(Acoeffs + Alen - 1, Acoeffs + Alen - 1, Acoeffs + Alen);
        Alen -= fmpz_is_zero(Acoeffs + Alen - 1);
    }

    /* from the fmpz_add: at most two junk coeffs past length */
    for (i = Alen; i < Alen + 2 && i < A->zpoly->alloc; i++)
        _fmpz_demote(Acoeffs + i);

    _fmpz_mpoly_set_length(A->zpoly, Alen, ctx->zctx);

    fmpz_clear(k);
    fmpz_clear(degmk);

    TMP_END;

    if (success)
    {
        fmpz_set(fmpq_numref(A->content), fmpq_numref(B->content));
        success = fmpz_pow_cache_mulpow_fmpz(fmpq_denref(A->content),
                                      fmpq_denref(B->content), deg, den_cache);
        fmpq_canonicalise(A->content);
    }
    else
    {
        fmpq_zero(A->content);
    }

    if (need_sort)
    {
        fmpz_mpoly_sort_terms(A->zpoly, ctx->zctx);
        fmpz_mpoly_combine_like_terms(A->zpoly, ctx->zctx);
    }

    fmpq_mpoly_reduce(A, ctx);

    FLINT_ASSERT(fmpq_mpoly_is_canonical(A, ctx));

    return success;
}

int fmpq_mpoly_evaluate_one_fmpq(fmpq_mpoly_t A,
                           const fmpq_mpoly_t B, slong var, const fmpq_t val,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    flint_bitcnt_t height;
    fmpz_pow_cache_t num_cache, den_cache;
    int success;

    if (B->zpoly->length == 0)
    {
        fmpq_mpoly_zero(A, ctx);
        return 1;
    }

    if (A == B)
    {
        int success;
        fmpq_mpoly_t T;
        fmpq_mpoly_init(T, ctx);
        success = fmpq_mpoly_evaluate_one_fmpq(T, B, var, val, ctx);
        fmpq_mpoly_swap(A, T, ctx);
        fmpq_mpoly_clear(T, ctx);
        return success;
    }

    fmpz_pow_cache_init(num_cache, fmpq_numref(val));
    fmpz_pow_cache_init(den_cache, fmpq_denref(val));
    height = fmpq_height_bits(val);

    if (B->zpoly->bits <= FLINT_BITS)
    {
        ulong deg = fmpq_mpoly_degree_si(B, var, ctx);
        success = (!_fmpz_pow_ui_is_not_feasible(height, deg)) &&
                  _fmpq_mpoly_evaluate_one_fmpq_sp(A, B, var,
                                               num_cache, den_cache, deg, ctx);
    }
    else
    {
        fmpz_t deg;
        fmpz_init(deg);
        fmpq_mpoly_degree_fmpz(deg, B, var, ctx);
        success = (!_fmpz_pow_fmpz_is_not_feasible(height, deg)) &&
                  _fmpq_mpoly_evaluate_one_fmpq_mp(A, B, var,
                                               num_cache, den_cache, deg, ctx);
        fmpz_clear(deg);
    }

    fmpz_pow_cache_clear(num_cache);
    fmpz_pow_cache_clear(den_cache);

    return success;
}

