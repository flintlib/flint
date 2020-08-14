/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

void fq_zech_mpoly_sub_fq_zech(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                              const fq_zech_t c, const fq_zech_mpoly_ctx_t ctx)
{
    slong i, N;
    slong Blen = B->length;

    if (fq_zech_is_zero(c, ctx->fqctx))
    {
        fq_zech_mpoly_set(A, B, ctx);
        return;
    }

    if (Blen == 0)
    {
        fq_zech_mpoly_set_fq_zech(A, c, ctx);
        FLINT_ASSERT(A->length > 0);
        fq_zech_neg(A->coeffs + 0, A->coeffs + 0, ctx->fqctx);
        return;
    }

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    if (mpoly_monomial_is_zero(B->exps + (Blen - 1)*N, N))
    {
        if (A != B)
        {
            fq_zech_mpoly_fit_length(A, B->length, ctx);
            fq_zech_mpoly_fit_bits(A, B->bits, ctx);
            A->bits = B->bits;

            for (i = 0; i < Blen - 1; i++)
                fq_zech_set(A->coeffs + i, B->coeffs + i, ctx->fqctx);

            mpoly_copy_monomials(A->exps, B->exps, Blen, N);

            _fq_zech_mpoly_set_length(A, B->length, ctx);

        }

        fq_zech_sub(A->coeffs + Blen - 1, B->coeffs + Blen - 1, c, ctx->fqctx);

        if (fq_zech_is_zero(A->coeffs + Blen - 1, ctx->fqctx))
            _fq_zech_mpoly_set_length(A, Blen - 1, ctx);
    }
    else
    {
        fq_zech_mpoly_fit_length(A, Blen + 1, ctx);

        if (A != B)
        {
            fq_zech_mpoly_fit_bits(A, B->bits, ctx);
            A->bits = B->bits;

            for (i = 0; i < Blen; i++)
                fq_zech_set(A->coeffs + i, B->coeffs + i, ctx->fqctx);

            mpoly_copy_monomials(A->exps, B->exps, Blen, N);
        } 

        mpoly_monomial_zero(A->exps + N*Blen, N);

        fq_zech_neg(A->coeffs + Blen, c, ctx->fqctx);

        _fq_zech_mpoly_set_length(A, B->length + 1, ctx);
    }
}
