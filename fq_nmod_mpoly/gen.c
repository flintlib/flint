/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, slong var,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;

    bits = mpoly_gen_bits_required(var, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(A, 1, ctx);
    fq_nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    fq_nmod_one(A->coeffs + 0, ctx->fqctx);
    if (bits <= FLINT_BITS)
        mpoly_gen_monomial_sp(A->exps, var, bits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(A->exps, var, bits, ctx->minfo);

    A->length = WORD(1);
}
