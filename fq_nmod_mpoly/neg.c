/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void _fq_nmod_mpoly_neg(fq_nmod_struct * Acoeff, ulong * Aexp,
                  const fq_nmod_struct * Bcoeff, const ulong * Bexp, slong Blen,
                                            slong N, const fq_nmod_ctx_t fqctx)
{
    slong i;

    for (i = 0; i < Blen; i++)
        fq_nmod_neg(Acoeff + i, Bcoeff + i, fqctx);

    if (Aexp != Bexp)
    {
        mpoly_copy_monomials(Aexp, Bexp, Blen, N);
    }
}

void fq_nmod_mpoly_neg(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;

    fq_nmod_mpoly_fit_length(A, B->length, ctx);
    fq_nmod_mpoly_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    _fq_nmod_mpoly_neg(A->coeffs, A->exps, B->coeffs, B->exps, B->length,
                                                                N, ctx->fqctx);
    _fq_nmod_mpoly_set_length(A, B->length, ctx);
}
