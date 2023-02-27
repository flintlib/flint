/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

void fq_zech_mpoly_scalar_mul_fq_zech(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_t B,
    const fq_zech_t c,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i, N;

    if (fq_zech_is_zero(c, ctx->fqctx))
    {
        fq_zech_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        if (fq_zech_is_one(c, ctx->fqctx))
            return;
    }
    else
    {
        fq_zech_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
        A->length = B->length;

        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        mpoly_copy_monomials(A->exps, B->exps, B->length, N);

        if (fq_zech_is_one(c, ctx->fqctx))
        {
            for (i = 0; i < B->length; i++)
                fq_zech_set(A->coeffs + i, B->coeffs + i, ctx->fqctx);

            return;
        }
    }

    for (i = 0; i < B->length; i++)
        fq_zech_mul(A->coeffs + i, B->coeffs + i, c, ctx->fqctx);
}
