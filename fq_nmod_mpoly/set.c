/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void _fq_nmod_mpoly_set(fq_nmod_struct * coeff1, ulong * exps1,
                const fq_nmod_struct * coeff2, const ulong * exps2, slong len2,
                                            slong N, const fq_nmod_ctx_t fqctx)
{
    slong i;

    if (coeff1 != coeff2)
    {
        for (i = 0; i < len2; i++)
            fq_nmod_set(coeff1 + i, coeff2 + i, fqctx);
    }

    if (exps1 != exps2)
    {
        mpoly_copy_monomials(exps1, exps2, len2, N);
    }
}

void fq_nmod_mpoly_set(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(A, B->length, ctx);
    fq_nmod_mpoly_fit_bits(A, B->bits, ctx);
    A->bits = B->bits;

    _fq_nmod_mpoly_set(A->coeffs, A->exps,
                                 B->coeffs, B->exps, B->length, N, ctx->fqctx);
    A->length = B->length;
}
