/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


static slong _fmpz_mpoly_scalar_fmma1(
    fmpz * Acoeffs, ulong * Aexps,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    const fmpz_t c,
    const fmpz * Dcoeffs, const ulong * Dexps, slong Dlen,
    const fmpz_t e,
    ulong maskhi)
{
    slong i = 0, j = 0, k = 0;

    while (i < Blen && j < Dlen)
    {
        if ((Bexps[i]^maskhi) > (Dexps[j]^maskhi))
        {
            fmpz_mul(Acoeffs + k, Bcoeffs + i, c);
            Aexps[k] = Bexps[i];
            i++;
            k++;
        }
        else if ((Bexps[i]^maskhi) == (Dexps[j]^maskhi))
        {
            fmpz_fmma(Acoeffs + k, Bcoeffs + i, c, Dcoeffs + j, e);
            Aexps[k] = Bexps[i];
            k += !fmpz_is_zero(Acoeffs + k);
            i++;
            j++;
        }
        else
        {
            fmpz_mul(Acoeffs + k, Dcoeffs + j, e);
            Aexps[k] = Dexps[j];
            j++;
            k++;
        }
    }

    while (i < Blen)
    {
        fmpz_mul(Acoeffs + k, Bcoeffs + i, c);
        Aexps[k] = Bexps[i];
        i++;
        k++;
    }

    while (j < Dlen)
    {
        fmpz_mul(Acoeffs + k, Dcoeffs + j, e);
        Aexps[k] = Dexps[j];
        j++;
        k++;
    }

    return k;
}


static slong _fmpz_mpoly_scalar_fmma(
    fmpz * Acoeffs, ulong * Aexps,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    const fmpz_t c,
    const fmpz * Dcoeffs, const ulong * Dexps, slong Dlen,
    const fmpz_t e,
    slong N,
    const ulong * cmpmask)
{
    slong i = 0, j = 0, k = 0;

    if (N == 1)
    {
        return _fmpz_mpoly_scalar_fmma1(Acoeffs, Aexps,
                                        Bcoeffs, Bexps, Blen, c,
                                        Dcoeffs, Dexps, Dlen, e, cmpmask[0]);
    }

    while (i < Blen && j < Dlen)
    {
        int cmp = mpoly_monomial_cmp(Bexps + i*N, Dexps + j*N, N, cmpmask);
        if (cmp > 0)
        {
            fmpz_mul(Acoeffs + k, Bcoeffs + i, c);
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            i++;
            k++;
        }
        else if (cmp == 0)
        {
            fmpz_fmma(Acoeffs + k, Bcoeffs + i, c, Dcoeffs + j, e);
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            k += !fmpz_is_zero(Acoeffs + k);
            i++;
            j++;
        }
        else
        {
            fmpz_mul(Acoeffs + k, Dcoeffs + j, e);
            mpoly_monomial_set(Aexps + k*N, Dexps + j*N, N);
            j++;
            k++;
        }
    }

    while (i < Blen)
    {
        fmpz_mul(Acoeffs + k, Bcoeffs + i, c);
        mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
        i++;
        k++;
    }

    while (j < Dlen)
    {
        fmpz_mul(Acoeffs + k, Dcoeffs + j, e);
        mpoly_monomial_set(Aexps + k*N, Dexps + j*N, N);
        j++;
        k++;
    }

   return k;
}

/* A = A*a + B*b */
void fmpz_mpoly_scalar_fmma_inplace(
    fmpz_mpoly_t A,
    const fmpz_t a,
    const fmpz_mpoly_t B,
    const fmpz_t b,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, s, new_len, N;
    slong Alen = A->length;
    slong Blen = B->length;
    ulong * Bexps, * cmpmask;
    int cmp, freeBexps;
    flint_bitcnt_t Abits;
    fmpz_mpoly_t T;
    TMP_INIT;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(Blen > 0);
    FLINT_ASSERT(!fmpz_is_zero(a));
    FLINT_ASSERT(!fmpz_is_zero(b));

    TMP_START;

    if (A->bits <= B->bits)
    {
        Abits = B->bits;
        if (A->bits < B->bits)
            fmpz_mpoly_repack_bits_inplace(A, Abits, ctx);
        N = mpoly_words_per_exp(Abits, ctx->minfo);
        Bexps = B->exps;
        freeBexps = 0;
    }
    else
    {
        Abits = A->bits;
        N = mpoly_words_per_exp(Abits, ctx->minfo);
        Bexps = (ulong *) flint_malloc(N*Blen*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits, Blen, ctx->minfo);
        freeBexps = 1;
    }

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    for (s = 0; s < Alen/4; s++)
    {
        cmp = mpoly_monomial_cmp(A->exps + N*(Alen - s - 1),
                                 Bexps + N*0, N, cmpmask);
        if (cmp >= 0)
        {
            s += (cmp == 0);
            goto doit;
        }
    }

    fmpz_mpoly_init3(T, Alen + Blen, Abits, ctx);
    T->length = _fmpz_mpoly_scalar_fmma(T->coeffs, T->exps,
                                  A->coeffs, A->exps, Alen, a,
                                  B->coeffs, Bexps, Blen, b, N, cmpmask);
    fmpz_mpoly_swap(A, T, ctx);
    fmpz_mpoly_clear(T, ctx);
    goto cleanup;

doit:

    FLINT_ASSERT(0 <= s && s <= Alen);

    FLINT_ASSERT(s == 0 || mpoly_monomial_cmp(A->exps + N*(Alen - s),
                                                Bexps + N*0, N, cmpmask) <= 0);

    FLINT_ASSERT(s == Alen || mpoly_monomial_cmp(A->exps + N*(Alen - s - 1),
                                                 Bexps + N*0, N, cmpmask) > 0);

    fmpz_mpoly_fit_length(A, Alen + Blen + s, ctx);
    mpoly_copy_monomials(A->exps + N*(Alen + Blen), A->exps + N*(Alen - s), s, N);
    _fmpz_vec_swap(A->coeffs + Alen + Blen, A->coeffs + Alen - s, s);

    if (!fmpz_is_one(a))
        _fmpz_vec_scalar_mul_fmpz(A->coeffs, A->coeffs, Alen - s, a);

    new_len = _fmpz_mpoly_scalar_fmma(
                  A->coeffs + Alen - s, A->exps + N*(Alen - s),
                  A->coeffs + (Alen + Blen), A->exps + N*(Alen + Blen), s, a,
                                        B->coeffs, Bexps, Blen, b, N, cmpmask);
    for (i = 0; i < s; i++)
        _fmpz_demote(A->coeffs + Alen + Blen + i);

    _fmpz_mpoly_set_length(A, Alen - s + new_len, ctx);

cleanup:

    if (freeBexps)
        flint_free(Bexps);

    TMP_END;

    return;
}

/* A = B*c + D*e */
void fmpz_mpoly_scalar_fmma(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_t c,
    const fmpz_mpoly_t D,
    const fmpz_t e,
    const fmpz_mpoly_ctx_t ctx)
{
    slong len, N;
    flint_bitcnt_t Abits;
    ulong * Bexps = B->exps, * Dexps = D->exps;
    ulong * cmpmask;
    int freeBexps = 0, freeDexps = 0;
    TMP_INIT;

    if (fmpz_mpoly_is_zero(B, ctx) || fmpz_is_zero(c))
    {
        fmpz_mpoly_scalar_mul_fmpz(A, D, e, ctx);
        return;
    }
    else if (fmpz_mpoly_is_zero(D, ctx) || fmpz_is_zero(e))
    {
        fmpz_mpoly_scalar_mul_fmpz(A, B, c, ctx);
        return;
    }
    else if (A == B)
    {
        if (A == D)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_add(t, c, e);
            fmpz_mpoly_scalar_mul_fmpz(A, A, t, ctx);
            fmpz_clear(t);
        }
        else
        {
            fmpz_mpoly_scalar_fmma_inplace(A, c, D, e, ctx);
        }
        return;
    }
    else if (A == D)
    {
        fmpz_mpoly_scalar_fmma_inplace(A, e, B, c, ctx);
        return;
    }

    Abits = FLINT_MAX(B->bits, D->bits);
    N = mpoly_words_per_exp(Abits, ctx->minfo);

    TMP_START;
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    if (Abits != B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits, B->length, ctx->minfo);
    }

    if (Abits != D->bits)
    {
        freeDexps = 1;
        Dexps = (ulong *) flint_malloc(N*D->length*sizeof(ulong));
        mpoly_repack_monomials(Dexps, Abits, D->exps, D->bits, D->length, ctx->minfo);
    }

    fmpz_mpoly_fit_length_reset_bits(A, B->length + D->length, Abits, ctx);

    len = _fmpz_mpoly_scalar_fmma(A->coeffs, A->exps,
                                   B->coeffs, Bexps, B->length, c,
                                   D->coeffs, Dexps, D->length, e, N, cmpmask);

    _fmpz_mpoly_set_length(A, len, ctx);

    if (freeBexps)
        flint_free(Bexps);

    if (freeDexps)
        flint_free(Dexps);

    TMP_END;
}
