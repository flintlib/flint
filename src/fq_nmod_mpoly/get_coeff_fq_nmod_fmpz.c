/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(
    fq_nmod_t c,
    const fq_nmod_mpoly_t A,
    fmpz * const * exp,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong index;
    index = mpoly_monomial_index_pfmpz(A->exps, A->bits, A->length,
                                                              exp, ctx->minfo);
    if (index < 0)
    {
        fq_nmod_zero(c, ctx->fqctx);
    }
    else
    {
        slong d = fq_nmod_ctx_degree(ctx->fqctx);
        FLINT_ASSERT(index < A->length);
        n_fq_get_fq_nmod(c, A->coeffs + d*index, ctx->fqctx);
    }
}
