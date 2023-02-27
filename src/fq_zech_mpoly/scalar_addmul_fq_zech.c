/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

static slong _fq_zech_mpoly_scalar_addmul_fq_zech(
    fq_zech_struct * Acoeffs, ulong * Aexps,
    fq_zech_struct * Bcoeffs, const ulong * Bexps, slong Blen,
    fq_zech_struct * Ccoeffs, const ulong * Cexps, slong Clen,
    const fq_zech_t d,
    slong N,
    const ulong * cmpmask,
    const fq_zech_ctx_t fqctx)
{
    slong i = 0, j = 0, k = 0;
    fq_zech_t p;

    FLINT_ASSERT(!fq_zech_is_zero(d, fqctx));

    fq_zech_init(p, fqctx);

    while (i < Blen && j < Clen)
    {
        int cmp = mpoly_monomial_cmp(Bexps + i*N, Cexps + j*N, N, cmpmask);

        if (cmp > 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            fq_zech_set(Acoeffs + k, Bcoeffs + i, fqctx);
            i++;
            k++;
        }
        else if (cmp == 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            fq_zech_mul(p, Ccoeffs + j, d, fqctx);
            fq_zech_add(Acoeffs + k, Bcoeffs + i, p, fqctx);
            k += !fq_zech_is_zero(Acoeffs + k, fqctx);
            i++;
            j++;
        }
        else
        {
            mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
            fq_zech_mul(Acoeffs + k, Ccoeffs + j, d, fqctx);
            j++;
            k++;
        }
    }

    while (i < Blen)
    {
        mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
        fq_zech_set(Acoeffs + k, Bcoeffs + i, fqctx);
        i++;
        k++;
    }

    while (j < Clen)
    {
        mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
        fq_zech_mul(Acoeffs + k, Ccoeffs + j, d, fqctx);
        j++;
        k++;
    }

    fq_zech_clear(p, fqctx);

    return k;
}

void fq_zech_mpoly_scalar_addmul_fq_zech(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    const fq_zech_mpoly_t C,
    const fq_zech_t d,
    const fq_zech_mpoly_ctx_t ctx)
{
    flint_bitcnt_t Abits;
    slong N;
    ulong * Bexps = B->exps, * Cexps = C->exps;
    ulong * cmpmask;
    int freeBexps = 0, freeCexps = 0;
    TMP_INIT;

    if (fq_zech_mpoly_is_zero(B, ctx))
    {
        fq_zech_mpoly_scalar_mul_fq_zech(A, C, d, ctx);
        return;
    }
    else if (fq_zech_mpoly_is_zero(C, ctx) || fq_zech_is_zero(d, ctx->fqctx))
    {
        fq_zech_mpoly_set(A, B, ctx);
        return;
    }

    TMP_START;
    Abits = FLINT_MAX(B->bits, C->bits);
    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    if (Abits != B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits, B->length, ctx->minfo);
    }

    if (Abits != C->bits)
    {
        freeCexps = 1;
        Cexps = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits, C->length, ctx->minfo);
    }

    if (A == B || A == C)
    {
        fq_zech_mpoly_t T;
        fq_zech_mpoly_init3(T, B->length + C->length, Abits, ctx);
        T->length = _fq_zech_mpoly_scalar_addmul_fq_zech(T->coeffs, T->exps, 
                                                  B->coeffs, Bexps, B->length,
                                                  C->coeffs, Cexps, C->length,
                                                    d, N, cmpmask, ctx->fqctx);
        fq_zech_mpoly_swap(A, T, ctx);
        fq_zech_mpoly_clear(T, ctx);
    }
    else
    {
        fq_zech_mpoly_fit_length_reset_bits(A, B->length + C->length, Abits, ctx);
        A->length = _fq_zech_mpoly_scalar_addmul_fq_zech(A->coeffs, A->exps, 
                                                  B->coeffs, Bexps, B->length,
                                                  C->coeffs, Cexps, C->length,
                                                    d, N, cmpmask, ctx->fqctx);
    }
      
    if (freeBexps)
        flint_free(Bexps);

    if (freeCexps)
        flint_free(Cexps);

    TMP_END;
}
