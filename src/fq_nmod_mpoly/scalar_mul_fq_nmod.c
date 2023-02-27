/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


void fq_nmod_mpoly_scalar_mul_n_fq(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const mp_limb_t * c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i;
    mp_limb_t * t;
    TMP_INIT;

    if (_n_fq_is_zero(c, d))
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        if (_n_fq_is_one(c, d))
            return;
    }
    else
    {
        slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

        fq_nmod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
        A->length = B->length;

        mpoly_copy_monomials(A->exps, B->exps, B->length, N);
        if (_n_fq_is_one(c, d))
        {
            _nmod_vec_set(A->coeffs, B->coeffs, d*B->length);
            return;
        }
    }

    TMP_START;

    t = (mp_limb_t *) TMP_ALLOC(d*N_FQ_MUL_ITCH*sizeof(mp_limb_t));

    for (i = 0; i < B->length; i++)
        _n_fq_mul(A->coeffs + d*i, B->coeffs + d*i, c, ctx->fqctx, t);

    TMP_END;
}


void fq_nmod_mpoly_scalar_mul_fq_nmod(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_t c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i;
    mp_limb_t * t;
    TMP_INIT;

    if (fq_nmod_is_zero(c, ctx->fqctx))
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        if (fq_nmod_is_one(c, ctx->fqctx))
            return;
    }
    else
    {
        slong N;

        fq_nmod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
        A->length = B->length;

        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        mpoly_copy_monomials(A->exps, B->exps, B->length, N);
        if (fq_nmod_is_one(c, ctx->fqctx))
        {
            _nmod_vec_set(A->coeffs, B->coeffs, d*B->length);
            return;
        }
    }

    TMP_START;

    t = (mp_limb_t *) TMP_ALLOC(d*(1 + N_FQ_MUL_ITCH)*sizeof(mp_limb_t));
    n_fq_set_fq_nmod(t, c, ctx->fqctx);

    for (i = 0; i < B->length; i++)
        _n_fq_mul(A->coeffs + d*i, B->coeffs + d*i, t, ctx->fqctx, t + d);

    TMP_END;
}
