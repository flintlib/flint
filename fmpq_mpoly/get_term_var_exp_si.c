/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

slong fmpq_mpoly_get_term_var_exp_si(const fmpq_mpoly_t A, slong i,
                                         slong var, const fmpq_mpoly_ctx_t ctx)
{
    slong N;

    if ((ulong) i >= (ulong) A->zpoly->length)
    {
        flint_throw(FLINT_ERROR,
                       "Index out of range in fmpq_mpoly_get_term_var_exp_si");
    }

    N = mpoly_words_per_exp(A->zpoly->bits, ctx->zctx->minfo);
    return mpoly_get_monomial_var_exp_si(A->zpoly->exps + N*i, var,
                                             A->zpoly->bits, ctx->zctx->minfo);
}
