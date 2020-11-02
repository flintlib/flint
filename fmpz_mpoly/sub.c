/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


slong _fmpz_mpoly_sub1(
    fmpz * Acoeffs, ulong * Aexps,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    const fmpz * Ccoeffs, const ulong * Cexps, slong Clen,
    ulong maskhi)
{
    slong i = 0, j = 0, k = 0;

    while (i < Blen && j < Clen)
    {
        if ((Bexps[i]^maskhi) > (Cexps[j]^maskhi))
        {
            Aexps[k] = Bexps[i];
            fmpz_set(Acoeffs + k, Bcoeffs + i);
            i++;
            k++;
        }
        else if ((Bexps[i]^maskhi) == (Cexps[j]^maskhi))
        {
            Aexps[k] = Bexps[i];
            fmpz_sub(Acoeffs + k, Bcoeffs + i, Ccoeffs + j);
            k += !fmpz_is_zero(Acoeffs + k);
            i++;
            j++;
        }
        else
        {
            Aexps[k] = Cexps[j];
            fmpz_neg(Acoeffs + k, Ccoeffs + j);
            j++;
            k++;
        }
    }

    while (i < Blen)
    {
        Aexps[k] = Bexps[i];
        fmpz_set(Acoeffs + k, Bcoeffs + i);
        i++;
        k++;
    }

    while (j < Clen)
    {
        Aexps[k] = Cexps[j];
        fmpz_neg(Acoeffs + k, Ccoeffs + j);
        j++;
        k++;
    }

   return k;
}

slong _fmpz_mpoly_sub(
    fmpz * Acoeffs, ulong * Aexps,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    const fmpz * Ccoeffs, const ulong * Cexps, slong Clen,
    slong N,
    const ulong * cmpmask)
{
    slong i = 0, j = 0, k = 0;

    if (N == 1)
    {
        return _fmpz_mpoly_sub1(Acoeffs, Aexps, Bcoeffs, Bexps, Blen,
                                             Ccoeffs, Cexps, Clen, cmpmask[0]);
    }

    while (i < Blen && j < Clen)
    {
        int cmp = mpoly_monomial_cmp(Bexps + i*N, Cexps + j*N, N, cmpmask);

        if (cmp > 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            fmpz_set(Acoeffs + k, Bcoeffs + i);
            i++;
            k++;
        }
        else if (cmp == 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            fmpz_sub(Acoeffs + k, Bcoeffs + i, Ccoeffs + j);
            k += !fmpz_is_zero(Acoeffs + k);
            i++;
            j++;
        }
        else
        {
            mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
            fmpz_neg(Acoeffs + k, Ccoeffs + j);
            j++;
            k++;
        }
    }

    while (i < Blen)
    {
        mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
        fmpz_set(Acoeffs + k, Bcoeffs + i);
        i++;
        k++;
    }

    while (j < Clen)
    {
        mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
        fmpz_neg(Acoeffs + k, Ccoeffs + j);
        j++;
        k++;
    }

    return k;
}

void fmpz_mpoly_sub_inplace(fmpz_mpoly_t A, const fmpz_mpoly_t B,
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
    T->length = _fmpz_mpoly_sub(T->coeffs, T->exps,
                                  A->coeffs, A->exps, Alen,
                                  B->coeffs, Bexps, Blen, N, cmpmask);
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

    new_len = _fmpz_mpoly_sub(A->coeffs + Alen - s, A->exps + N*(Alen - s),
                     A->coeffs + (Alen + Blen), A->exps + N*(Alen + Blen), s,
                                           B->coeffs, Bexps, Blen, N, cmpmask);
    for (i = 0; i < s; i++)
        _fmpz_demote(A->coeffs + Alen + Blen + i);

    _fmpz_mpoly_set_length(A, Alen - s + new_len, ctx);

cleanup:

    if (freeBexps)
        flint_free(Bexps);

    TMP_END;

    return;
}

void fmpz_mpoly_sub(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                              const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
{
    slong Alen, N;
    flint_bitcnt_t Abits;
    ulong * Bexps = B->exps, * Cexps = C->exps;
    ulong * cmpmask;
    int freeBexps = 0, freeCexps = 0;
    TMP_INIT;

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        fmpz_mpoly_neg(A, C, ctx);
        return;
    }
    else if (fmpz_mpoly_is_zero(C, ctx))
    {
        fmpz_mpoly_set(A, B, ctx);
        return;
    }
    else if (A == B)
    {
        if (A == C)
            fmpz_mpoly_zero(A, ctx);
        else
            fmpz_mpoly_sub_inplace(A, C, ctx);
        return;
    }
    else if (A == C)
    {
        fmpz_mpoly_sub_inplace(A, B, ctx);
        _fmpz_vec_neg(A->coeffs, A->coeffs, A->length);
        return;
    }

    TMP_START;
    Abits = FLINT_MAX(B->bits, C->bits);
    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    if (Abits > B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits, B->length, ctx->minfo);
    }

    if (Abits > C->bits)
    {
        freeCexps = 1;
        Cexps = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits, C->length, ctx->minfo);
    }

    fmpz_mpoly_fit_length_reset_bits(A, B->length + C->length, Abits, ctx);

    Alen = _fmpz_mpoly_sub(A->coeffs, A->exps,
                                  B->coeffs, Bexps, B->length,
                                  C->coeffs, Cexps, C->length, N, cmpmask);
    _fmpz_mpoly_set_length(A, Alen, ctx);

    if (freeBexps)
        flint_free(Bexps);

    if (freeCexps)
        flint_free(Cexps);

    TMP_END;
}
