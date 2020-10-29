/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_is_gen(const fq_nmod_mpoly_t A,
                                      slong var, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);

    if (A->length != 1)
        return 0;

    if (!_n_fq_is_one(A->coeffs + d*0, d))
        return 0;

    return mpoly_is_gen(A->exps, var, A->bits, ctx->minfo);
}
