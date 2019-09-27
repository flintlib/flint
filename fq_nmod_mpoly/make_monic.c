/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_make_monic(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_t c;

    if (B->length == 0)
    {
        flint_throw(FLINT_ERROR, "Zero polynomial in fq_nmod_mpoly_make_monic");
    }

    fq_nmod_init(c, ctx->fqctx);
    fq_nmod_inv(c, B->coeffs + 0, ctx->fqctx);
    fq_nmod_mpoly_scalar_mul_fq_nmod(A, B, c, ctx);
    fq_nmod_clear(c, ctx->fqctx);
}
