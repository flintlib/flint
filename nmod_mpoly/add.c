/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

slong _nmod_mpoly_add1(
    mp_limb_t * Acoeffs, ulong * Aexps,
    const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    const mp_limb_t * Ccoeffs, const ulong * Cexps, slong Clen,
    ulong maskhi,
    nmod_t fctx)
{
    slong i = 0, j = 0, k = 0;

    while (i < Blen && j < Clen)
    {
        if ((Bexps[i]^maskhi) > (Cexps[j]^maskhi))
        {
            Aexps[k] = Bexps[i];
            Acoeffs[k] = Bcoeffs[i];
            i++;
        }
        else if ((Bexps[i]^maskhi) == (Cexps[j]^maskhi))
        {
            Aexps[k] = Bexps[i];
            Acoeffs[k] = nmod_add(Bcoeffs[i], Ccoeffs[j], fctx);
            k -= (Acoeffs[k] == 0);
            i++;
            j++;
        }
        else
        {
            Acoeffs[k] = Ccoeffs[j];
            Aexps[k] = Cexps[j];
            j++;         
        }
        k++;
    }

    while (i < Blen)
    {
        Aexps[k] = Bexps[i];
        Acoeffs[k] = Bcoeffs[i];
        i++;
        k++;
    }

    while (j < Clen)
    {
        Aexps[k] = Cexps[j];
        Acoeffs[k] = Ccoeffs[j];
        j++;
        k++;
    }

    return k;
}

slong _nmod_mpoly_add(mp_limb_t * Acoeffs,       ulong * Aexps,
                const mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
                const mp_limb_t * Ccoeffs, const ulong * Cexps, slong Clen,
                   slong N, const ulong * cmpmask, nmod_t fctx)
{
    slong i = 0, j = 0, k = 0;

    if (N == 1)
        return _nmod_mpoly_add1(Acoeffs, Aexps, Bcoeffs, Bexps, Blen,
                                       Ccoeffs, Cexps, Clen, cmpmask[0], fctx);

    while (i < Blen && j < Clen)
    {
        int cmp = mpoly_monomial_cmp(Bexps + i*N, Cexps + j*N, N, cmpmask);

        if (cmp > 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            Acoeffs[k] = Bcoeffs[i];
            i++;
        }
        else if (cmp == 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            Acoeffs[k] = nmod_add(Bcoeffs[i], Ccoeffs[j], fctx);
            k -= (Acoeffs[k] == 0);
            i++;
            j++;
        }
        else
        {
            Acoeffs[k] = Ccoeffs[j];
            mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
            j++;         
        }
        k++;
    }

    while (i < Blen)
    {
        mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
        Acoeffs[k] = Bcoeffs[i];
        i++;
        k++;
    }

    while (j < Clen)
    {
        mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
        Acoeffs[k] = Ccoeffs[j];
        j++;
        k++;
    }

    return k;
}

void nmod_mpoly_add(nmod_mpoly_t A, const nmod_mpoly_t B,
                              const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
{
    slong Abits, N;
    ulong * Bexps = B->exps, * Cexps = C->exps;
    ulong * cmpmask;
    int freeBexps = 0, freeCexps = 0;
    TMP_INIT;

    if (B->length == 0)
    {
        nmod_mpoly_set(A, C, ctx);
        return;
    }
    else if (C->length == 0)
    {
        nmod_mpoly_set(A, B, ctx);
        return;
    }

    TMP_START;
    Abits = FLINT_MAX(B->bits, C->bits);
    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    if (Abits != B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits,
                                                    B->length, ctx->minfo);
    }

    if (Abits != C->bits)
    {
        freeCexps = 1;
        Cexps = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits,
                                                    C->length, ctx->minfo);
    }

    if (A == B || A == C)
    {
        nmod_mpoly_t T;
        nmod_mpoly_init3(T, B->length + C->length, Abits, ctx);
        T->length = _nmod_mpoly_add(T->coeffs, T->exps, 
                                    B->coeffs, Bexps, B->length,
                                    C->coeffs, Cexps, C->length,
                                                      N, cmpmask, ctx->mod);
        nmod_mpoly_swap(A, T, ctx);
        nmod_mpoly_clear(T, ctx);
    }
    else
    {
        nmod_mpoly_fit_length_reset_bits(A, B->length + C->length, Abits, ctx);
        A->length = _nmod_mpoly_add(A->coeffs, A->exps, 
                                    B->coeffs, Bexps, B->length,
                                    C->coeffs, Cexps, C->length,
                                                      N, cmpmask, ctx->mod);
    }

    if (freeBexps)
        flint_free(Bexps);

    if (freeCexps)
        flint_free(Cexps);

    TMP_END;
}
