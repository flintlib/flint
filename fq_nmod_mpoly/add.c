/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

slong _fq_nmod_mpoly_add(
    mp_limb_t * Acoeffs, ulong * Aexps,
    mp_limb_t * Bcoeffs, const ulong * Bexps, slong Blen,
    mp_limb_t * Ccoeffs, const ulong * Cexps, slong Clen,
    slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx)
{
    slong d = fq_nmod_ctx_degree(fqctx);
    nmod_t mod = fq_nmod_ctx_mod(fqctx);
    slong i = 0, j = 0, k = 0;

    while (i < Blen && j < Clen)
    {
        int cmp = mpoly_monomial_cmp(Bexps + N*i, Cexps + N*j, N, cmpmask);

        if (cmp > 0)
        {
            mpoly_monomial_set(Aexps + N*k, Bexps + N*i, N);
            _n_fq_set(Acoeffs + d*k, Bcoeffs + d*i, d);
            i++;
        }
        else if (cmp == 0)
        {
            mpoly_monomial_set(Aexps + N*k, Bexps + N*i, N);
            _n_fq_add(Acoeffs + d*k, Bcoeffs + d*i, Ccoeffs + d*j, d, mod);
            k -= _n_fq_is_zero(Acoeffs + d*k, d);
            i++;
            j++;
        }
        else
        {
            mpoly_monomial_set(Aexps + N*k, Cexps + N*j, N);
            _n_fq_set(Acoeffs + d*k, Ccoeffs + d*j, d);
            j++;
        }
        k++;
    }

    while (i < Blen)
    {
        mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
        _n_fq_set(Acoeffs + d*k, Bcoeffs + d*i, d);
        i++;
        k++;
    }

    while (j < Clen)
    {
        mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
        _n_fq_set(Acoeffs + d*k, Ccoeffs + d*j, d);
        j++;
        k++;
    }

    return k;
}

void fq_nmod_mpoly_add(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_t C,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong Alen;
    ulong * Bexps = B->exps, * Cexps = C->exps;
    flint_bitcnt_t Abits = FLINT_MAX(B->bits, C->bits);
    slong N = mpoly_words_per_exp(Abits, ctx->minfo);
    ulong * cmpmask;
    int freeBexps = 0, freeCexps = 0;
    TMP_INIT;

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_set(A, C, ctx);
        return;
    }
    else if (fq_nmod_mpoly_is_zero(C, ctx))
    {
        fq_nmod_mpoly_set(A, B, ctx);
        return;
    }

    TMP_START;
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    if (Abits > B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits,
                                                    B->length, ctx->minfo);
    }

    if (Abits > C->bits)
    {
        freeCexps = 1;
        Cexps = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits,
                                                    C->length, ctx->minfo);
    }

    if (A == B || A == C)
    {
        fq_nmod_mpoly_t T;
        fq_nmod_mpoly_init(T, ctx);
        fq_nmod_mpoly_fit_length_reset_bits(T, B->length + C->length, Abits, ctx);
        Alen = _fq_nmod_mpoly_add(T->coeffs, T->exps, 
                                      B->coeffs, Bexps, B->length,
                                      C->coeffs, Cexps, C->length,
                                                       N, cmpmask, ctx->fqctx);
        fq_nmod_mpoly_swap(T, A, ctx);
        fq_nmod_mpoly_clear(T, ctx);
    }
    else
    {
        fq_nmod_mpoly_fit_length_reset_bits(A, B->length + C->length, Abits, ctx);
        Alen = _fq_nmod_mpoly_add(A->coeffs, A->exps, 
                                      B->coeffs, Bexps, B->length,
                                      C->coeffs, Cexps, C->length,
                                                       N, cmpmask, ctx->fqctx);
    }

    _fq_nmod_mpoly_set_length(A, Alen, ctx);
      
    if (freeBexps)
        flint_free(Bexps);

    if (freeCexps)
        flint_free(Cexps);

    TMP_END;
}
