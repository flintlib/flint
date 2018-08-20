/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void _fq_nmod_mpoly_neg(fq_nmod_struct * coeff1, ulong * exp1,
                  const fq_nmod_struct * coeff2, const ulong * exp2, slong len,
                                            slong N, const fq_nmod_ctx_t fqctx)
{
    slong i;

    for (i = 0; i < len; i++)
        fq_nmod_neg(coeff1 + i, coeff2 + i, fqctx);

    if (exp1 != exp2)
    {
        for (i = 0; i < len; i++)
            mpoly_monomial_set(exp1 + N*i, exp2 + N*i, N);
    }
}

void fq_nmod_mpoly_neg(fq_nmod_mpoly_t poly1, const fq_nmod_mpoly_t poly2,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    N = mpoly_words_per_exp(poly2->bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(poly1, poly2->length, ctx);
    fq_nmod_mpoly_fit_bits(poly1, poly2->bits, ctx);

    _fq_nmod_mpoly_neg(poly1->coeffs, poly1->exps,
                   poly2->coeffs, poly2->exps, poly2->length, N, ctx->fqctx);

    _fq_nmod_mpoly_set_length(poly1, poly2->length, ctx);
    poly1->bits = poly2->bits;
}
