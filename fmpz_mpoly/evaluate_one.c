/*
    Copyright (C) 2018, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_pow_cache_init(fmpz_pow_cache_t T, const fmpz_t val)
{
    fmpz_init(T->tmp);
    T->alloc = 10;
    T->powers = _fmpz_vec_init(T->alloc);
    fmpz_one(T->powers + 0);
    fmpz_set(T->powers + 1, val);
    T->length = 2;
}

void fmpz_pow_cache_clear(fmpz_pow_cache_t T)
{
    fmpz_clear(T->tmp);
    _fmpz_vec_clear(T->powers, T->alloc);
}

/* a = b * val^k */
int fmpz_pow_cache_mulpow_ui(
    fmpz_t a,
    const fmpz_t b,
    ulong k,
    fmpz_pow_cache_t T)
{
    slong i;

    if (k > 100)
    {
        fmpz_pow_ui(T->tmp, T->powers + 1, k);
        fmpz_mul(a, b, T->tmp);
        return 1;
    }

    if (k >= T->length)
    {
        if (k + 1 >= T->alloc)
        {
            slong new_alloc = FLINT_MAX(k + 1, 2*T->alloc);
            T->powers = FLINT_ARRAY_REALLOC(T->powers, new_alloc, fmpz);
            for (i = T->alloc; i < new_alloc; i++)
                fmpz_init(T->powers + i);

            T->alloc = new_alloc;
        }

        do {
            fmpz_mul(T->powers + T->length, T->powers + T->length - 1,
                                                                T->powers + 1);
            T->length++;
        } while (k >= T->length);
    }

    fmpz_mul(a, b, T->powers + k);
    return 1;
}

int fmpz_pow_cache_mulpow_fmpz(
    fmpz_t a,
    const fmpz_t b,
    const fmpz_t k,
    fmpz_pow_cache_t T)
{
    if (fmpz_abs_fits_ui(k))
        return fmpz_pow_cache_mulpow_ui(a, b, fmpz_get_ui(k), T);

    if (!fmpz_pow_fmpz(T->tmp, T->powers + 1, k))
        return 0;

    fmpz_mul(a, b, T->tmp);
    return 1;
}

/* exponents of B are not multiprecision */
static int _fmpz_mpoly_evaluate_one_fmpz_sp(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong var,
    fmpz_pow_cache_t cache,
    const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i, N, off, shift;
    ulong * cmpmask, * one;
    slong Blen = B->length;
    const fmpz * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t bits = B->bits;
    slong Alen;
    fmpz * Acoeffs;
    ulong * Aexps;
    ulong mask, k;
    int need_sort = 0, cmp;
    TMP_INIT;

    FLINT_ASSERT(B->bits <= FLINT_BITS);

    TMP_START;

    if (A != B)
        fmpz_mpoly_fit_length_reset_bits(A, Blen, bits, ctx);

    Acoeffs = A->coeffs;
    Aexps = A->exps;

    mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    cmpmask = TMP_ARRAY_ALLOC(2*N, ulong);
    one = cmpmask + N;
    mpoly_gen_monomial_offset_shift_sp(one, &off, &shift, var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        k = (Bexps[N*i + off] >> shift) & mask;
        success = fmpz_pow_cache_mulpow_ui(Acoeffs + Alen, Bcoeffs + i, k, cache);
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
    for (i = Alen; i < Alen + 2 && i < A->alloc; i++)
        _fmpz_demote(Acoeffs + i);

    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;

    if (need_sort)
    {
        fmpz_mpoly_sort_terms(A, ctx);
        fmpz_mpoly_combine_like_terms(A, ctx);
    }

    FLINT_ASSERT(fmpz_mpoly_is_canonical(A, ctx));

    return success;
}

/* exponents of B are multiprecision */
static int _fmpz_mpoly_evaluate_one_fmpz_mp(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong var,
    fmpz_pow_cache_t cache,
    const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i, N, off;
    ulong * cmpmask, * one, * tmp;
    slong Blen = B->length;
    const fmpz * Bcoeffs = B->coeffs;
    const ulong * Bexps = B->exps;
    flint_bitcnt_t bits = B->bits;
    slong Alen;
    fmpz * Acoeffs;
    ulong * Aexps;
    fmpz_t k;
    int need_sort = 0, cmp;
    TMP_INIT;

    FLINT_ASSERT((B->bits % FLINT_BITS) == 0);

    TMP_START;

    fmpz_init(k);
    
    if (A != B)
        fmpz_mpoly_fit_length_reset_bits(A, Blen, bits, ctx);

    Acoeffs = A->coeffs;
    Aexps = A->exps;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong *) TMP_ALLOC(3*N*sizeof(ulong));
    cmpmask = one + N;
    tmp = cmpmask + N;
    off = mpoly_gen_monomial_offset_mp(one, var, bits, ctx->minfo);
    mpoly_get_cmpmask(cmpmask, N, bits, ctx->minfo);

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        fmpz_set_ui_array(k, Bexps + N*i + off, bits/FLINT_BITS);
        success = fmpz_pow_cache_mulpow_fmpz(Acoeffs + Alen, Bcoeffs + i, k, cache);
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
    for (i = Alen; i < Alen + 2 && i < A->alloc; i++)
        _fmpz_demote(Acoeffs + i);

    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_clear(k);

    TMP_END;

    if (need_sort)
    {
        fmpz_mpoly_sort_terms(A, ctx);
        fmpz_mpoly_combine_like_terms(A, ctx);
    }

    FLINT_ASSERT(fmpz_mpoly_is_canonical(A, ctx));

    return success;
}

int fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                       slong var, const fmpz_t val, const fmpz_mpoly_ctx_t ctx)
{
    int success;
    fmpz_pow_cache_t T;

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        fmpz_mpoly_zero(A, ctx);
        return 1;
    }

    fmpz_pow_cache_init(T, val);

    if (B->bits <= FLINT_BITS)
        success = _fmpz_mpoly_evaluate_one_fmpz_sp(A, B, var, T, ctx);
    else
        success = _fmpz_mpoly_evaluate_one_fmpz_mp(A, B, var, T, ctx);

    fmpz_pow_cache_clear(T);

    return success;
}
