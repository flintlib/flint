/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/* setform copies the exponents and zeros the coefficients */
void nmod_mpoly_setform(nmod_mpoly_t A, nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = B->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    nmod_mpoly_fit_length_reset_bits(A, B->length, bits, ctx);
    mpoly_copy_monomials(A->exps, B->exps, B->length, N);
    _nmod_vec_zero(A->coeffs, B->length);
    A->length = B->length;
}

void nmod_mpolyu_setform(nmod_mpolyu_t A, nmod_mpolyu_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        nmod_mpoly_setform(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}

void nmod_mpoly_setform_mpolyn(nmod_mpoly_t A, nmod_mpolyn_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong N;

    FLINT_ASSERT(A->bits == B->bits);

    nmod_mpoly_fit_length(A, B->length, ctx);
    N = mpoly_words_per_exp(B->bits, ctx->minfo);
    for (i = 0; i < B->length; i++)
    {
        A->coeffs[i] = UWORD(0);
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
    }
    A->length = B->length;
}

void nmod_mpolyu_setform_mpolyun(nmod_mpolyu_t A, nmod_mpolyun_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i;

    nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        FLINT_ASSERT((B->coeffs + i)->bits == B->bits);
        nmod_mpoly_setform_mpolyn(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}
