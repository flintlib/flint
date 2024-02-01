/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_term_exp_fits_ui(const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_ui(A->exps, A->bits, i, ctx->minfo);
}

int fq_nmod_mpoly_term_exp_fits_si(const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
{
    return A->bits <= FLINT_BITS ? 1
                     : mpoly_term_exp_fits_si(A->exps, A->bits, i, ctx->minfo);
}
