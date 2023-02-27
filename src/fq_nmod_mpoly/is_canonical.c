/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_is_canonical(
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong i;

    if (d*A->length > A->coeffs_alloc)
        return 0;

    if (N*A->length > A->exps_alloc)
        return 0;

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!n_fq_is_canonical(A->coeffs + d*i, ctx->fqctx))
            return 0;

        if (_n_fq_is_zero(A->coeffs + d*i, d))
            return 0;
    }

    return 1;
}

void fq_nmod_mpoly_assert_canonical(
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong i;

    if (d*A->length > A->coeffs_alloc)
        flint_throw(FLINT_ERROR, "Polynomial coefficient allocation is bad");

    if (N*A->length > A->exps_alloc)
        flint_throw(FLINT_ERROR, "Polynomial exponent allocation is bad");

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents out of order");

    for (i = 0; i < A->length; i++)
    {
        if (!n_fq_is_canonical(A->coeffs + d*i, ctx->fqctx))
            flint_throw(FLINT_ERROR, "Polynomial has a bad coefficient");

        if (_n_fq_is_zero(A->coeffs + d*i, d))
            flint_throw(FLINT_ERROR, "Polynomial has a zero coefficient");
    }
}
