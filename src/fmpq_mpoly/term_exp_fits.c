/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fmpq_mpoly.h"

int fmpq_mpoly_term_exp_fits_ui(const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
{
    return A->zpoly->bits <= FLINT_BITS ? 1
                               : mpoly_term_exp_fits_ui(A->zpoly->exps,
                                          A->zpoly->bits, i, ctx->zctx->minfo);
}

int fmpq_mpoly_term_exp_fits_si(const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
{
    return A->zpoly->bits <= FLINT_BITS ? 1
                               : mpoly_term_exp_fits_si(A->zpoly->exps,
                                          A->zpoly->bits, i, ctx->zctx->minfo);
}
