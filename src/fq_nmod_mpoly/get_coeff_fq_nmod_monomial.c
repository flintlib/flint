/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_get_coeff_fq_nmod_monomial(
    fq_nmod_t c,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t M,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong index, N;
    ulong * cmpmask, * pexp;
    int exists;
    TMP_INIT;

    if (M->length != WORD(1))
    {
        flint_throw(FLINT_ERROR,
                 "M not monomial in fq_nmod_mpoly_get_coeff_fq_nmod_monomial");
    }

    TMP_START;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    pexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));

    mpoly_get_cmpmask(cmpmask, N, A->bits, ctx->minfo);

    if (M->bits == A->bits)
    {
        mpoly_monomial_set(pexp, M->exps + N*0, N);
    }
    else
    {
        int could_repack = mpoly_repack_monomials(pexp, A->bits,
                                        M->exps + N*0, M->bits, 1, ctx->minfo);
        if (!could_repack)
        {
            FLINT_ASSERT(M->bits > A->bits);
            fq_nmod_zero(c, ctx->fqctx);
            goto clean_up;
        }
    }

    exists = mpoly_monomial_exists(&index, A->exps, pexp, A->length, N, cmpmask);

    if (!exists)
    {
        fq_nmod_zero(c, ctx->fqctx);
    }
    else
    {
        slong d = fq_nmod_ctx_degree(ctx->fqctx);
        FLINT_ASSERT(index < A->length);
        n_fq_get_fq_nmod(c, A->coeffs + d*index, ctx->fqctx);
    }

clean_up:

    TMP_END;
    return;
}
