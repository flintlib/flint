/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "n_poly.h"
#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_is_monic(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    return A->length > 0 &&
           _n_fq_is_one(A->coeffs + 0, fq_nmod_ctx_degree(ctx->fqctx));
}
