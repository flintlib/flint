/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

/* TODO move this and make it faster */
static int fq_zech_get_ui(ulong * res, const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    int success;
    nmod_poly_t asdf;

    nmod_poly_init_mod(asdf, fq_zech_ctx_modulus(ctx)->mod);

    fq_zech_get_nmod_poly(asdf, op, ctx);

    success = (asdf->length <= 1);
    if (asdf->length == 1)
        *res = asdf->coeffs[0];
    else
        *res = 0;

    nmod_poly_clear(asdf);
    return success;
}

int _fq_zech_mpoly_get_nmod_mpoly(
    nmod_mpoly_t s,
    const nmod_mpoly_ctx_t sctx,
    const fq_zech_mpoly_t t,
    const fq_zech_mpoly_ctx_t tctx)
{
    slong i, N;

    FLINT_ASSERT(sctx->minfo->ord == tctx->minfo->ord);
    FLINT_ASSERT(sctx->minfo->nvars == tctx->minfo->nvars);

    nmod_mpoly_fit_length_reset_bits(s, t->length, t->bits, sctx);
    s->length = t->length;

    N = mpoly_words_per_exp(t->bits, tctx->minfo);
    mpoly_copy_monomials(s->exps, t->exps, t->length, N);

    for (i = 0; i < t->length; i++)
    {
        if (!fq_zech_get_ui(s->coeffs + i, t->coeffs + i, tctx->fqctx))
            return 0;
    }

    return 1;
}

void _fq_zech_mpoly_set_nmod_mpoly(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t Actx,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t Bctx)
{
    slong i, N;

    FLINT_ASSERT(Actx->minfo->ord == Bctx->minfo->ord);
    FLINT_ASSERT(Actx->minfo->nvars == Bctx->minfo->nvars);

    fq_zech_mpoly_fit_length_reset_bits(A, B->length, B->bits, Actx);
    A->length = B->length;

    N = mpoly_words_per_exp(B->bits, Bctx->minfo);
    mpoly_copy_monomials(A->exps, B->exps, B->length, N);

    for (i = 0; i < B->length; i++)
        fq_zech_set_ui(A->coeffs + i, B->coeffs[i], Actx->fqctx);
}

