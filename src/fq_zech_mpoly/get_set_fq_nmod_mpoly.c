/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"
#include "fq_nmod_mpoly_factor.h"


void _fq_zech_mpoly_get_fq_nmod_mpoly(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctxA,
    const fq_zech_mpoly_t B,
    const fq_zech_mpoly_ctx_t ctxB)
{
    slong d = fq_nmod_ctx_degree(ctxA->fqctx);
    flint_bitcnt_t bits = B->bits;
    slong i, N = mpoly_words_per_exp(bits, ctxB->minfo);
    fq_nmod_t t;

    fq_nmod_init(t, ctxA->fqctx);

    FLINT_ASSERT(ctxA->minfo->ord == ctxB->minfo->ord);
    FLINT_ASSERT(ctxA->minfo->nvars == ctxB->minfo->nvars);

    fq_nmod_mpoly_fit_length_reset_bits(A, B->length, bits, ctxA);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
    {
        fq_zech_get_fq_nmod(t, B->coeffs + i, ctxB->fqctx);
        n_fq_set_fq_nmod(A->coeffs + d*i, t, ctxA->fqctx);
    }
    mpoly_copy_monomials(A->exps, B->exps, B->length, N);

    fq_nmod_clear(t, ctxA->fqctx);
}

void _fq_zech_mpoly_set_fq_nmod_mpoly(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctxA,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctxB)
{
    slong d = fq_nmod_ctx_degree(ctxB->fqctx);
    flint_bitcnt_t bits = B->bits;
    slong i, N = mpoly_words_per_exp(bits, ctxB->minfo);
    fq_nmod_t t;

    fq_nmod_init(t, ctxB->fqctx);

    FLINT_ASSERT(ctxA->minfo->ord == ctxB->minfo->ord);
    FLINT_ASSERT(ctxA->minfo->nvars == ctxB->minfo->nvars);

    fq_zech_mpoly_fit_length_reset_bits(A, B->length, bits, ctxA);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
    {
        n_fq_get_fq_nmod(t, B->coeffs + d*i, ctxB->fqctx);
        fq_zech_set_fq_nmod(A->coeffs + i, t, ctxA->fqctx);
    }
    mpoly_copy_monomials(A->exps, B->exps, B->length, N);

    fq_nmod_clear(t, ctxB->fqctx);
}
