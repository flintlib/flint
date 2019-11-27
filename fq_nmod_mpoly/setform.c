/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/* setform copies the exponents and zeros the coefficients */
void fq_nmod_mpoly_setform(fq_nmod_mpoly_t A, fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpoly_set(A, B, ctx);
    for (i = 0; i < A->length; i++)
    {
        fq_nmod_zero(A->coeffs + i, ctx->fqctx);
    }
}

void fq_nmod_mpolyu_setform(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    fq_nmod_mpolyu_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mpoly_setform(A->coeffs + i, B->coeffs + i, ctx);
        A->exps[i] = B->exps[i];
    }
    A->length = B->length;
}
