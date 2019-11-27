/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_get_term_monomial(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A,
                                        slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t bits = A->bits;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR,
                      "Index out of range in fq_nmod_mpoly_get_term_monomial");
    }

    fq_nmod_mpoly_fit_length(M, WORD(1), ctx);
    fq_nmod_mpoly_fit_bits(M, bits, ctx);
    M->bits = bits;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    mpoly_monomial_set(M->exps + N*0, A->exps + N*i, N);
    fq_nmod_one(M->coeffs + 0, ctx->fqctx);
    _fq_nmod_mpoly_set_length(M, 1, ctx);
}
