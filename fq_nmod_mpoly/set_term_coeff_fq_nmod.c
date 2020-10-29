/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_set_term_coeff_fq_nmod(fq_nmod_mpoly_t A,
                     slong i, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);

    if (i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "fq_nmod_mpoly_set_term_coeff_fq_nmod: index is out of range");
    }

    n_fq_set_fq_nmod(A->coeffs + d*i, c, ctx->fqctx);
}
