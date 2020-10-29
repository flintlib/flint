/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void _fq_nmod_mpoly_push_exp_ffmpz(
    fq_nmod_mpoly_t A,
    const fmpz * exp,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    slong old_length = A->length;
    flint_bitcnt_t exp_bits;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fq_nmod_mpoly_fit_length_fit_bits(A, old_length + 1, exp_bits, ctx);

    A->length = old_length + 1;
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_set_monomial_ffmpz(A->exps + N*old_length, exp, A->bits, ctx->minfo);
}


void _fq_nmod_mpoly_push_exp_pfmpz(
    fq_nmod_mpoly_t A,
    fmpz * const * exp,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;
    slong old_length = A->length;
    flint_bitcnt_t exp_bits;

    exp_bits = mpoly_exp_bits_required_pfmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fq_nmod_mpoly_fit_length_fit_bits(A, old_length + 1, exp_bits, ctx);

    A->length = old_length + 1;
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_set_monomial_pfmpz(A->exps + N*old_length, exp, A->bits, ctx->minfo);
}


void fq_nmod_mpoly_push_term_fq_nmod_fmpz(
    fq_nmod_mpoly_t A,
    const fq_nmod_t c,
    fmpz * const * exp,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d;
    _fq_nmod_mpoly_push_exp_pfmpz(A, exp, ctx);
    d = fq_nmod_ctx_degree(ctx->fqctx);
    FLINT_ASSERT(A->length > 0);
    n_fq_set_fq_nmod(A->coeffs + d*(A->length - 1), c, ctx->fqctx);
}
