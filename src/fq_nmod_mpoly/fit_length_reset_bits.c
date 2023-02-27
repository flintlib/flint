/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void fq_nmod_mpoly_fit_length_reset_bits(
    fq_nmod_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    _fq_nmod_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc, d,
                              &A->exps, &A->exps_alloc, N, len);
    A->bits = bits;
}
