/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "fq_nmod.h"
#include "mpoly.h"
#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_deflate(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fmpz * shift,
    const fmpz * stride,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong bits = B->bits;
    slong NA = mpoly_words_per_exp(bits, ctx->minfo);

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        slong new_alloc = NA*A->length;
        ulong * texps = flint_malloc(new_alloc*sizeof(ulong));
        mpoly_monomials_deflate(texps, bits, B->exps, B->bits, B->length,
                                                    shift, stride, ctx->minfo);
        flint_free(A->exps);
        A->exps = texps;
        A->bits = bits;
        A->exps_alloc = new_alloc;
    }
    else
    {
        fq_nmod_mpoly_fit_length_reset_bits(A, B->length, bits, ctx);
        _nmod_vec_set(A->coeffs, B->coeffs, d*B->length);
        mpoly_monomials_deflate(A->exps, bits, B->exps, B->bits, B->length,
                                                    shift, stride, ctx->minfo);
        _fq_nmod_mpoly_set_length(A, B->length, ctx);
    }

    if (ctx->minfo->ord != ORD_LEX)
        fq_nmod_mpoly_sort_terms(A, ctx);
}
