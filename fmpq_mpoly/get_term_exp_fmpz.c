/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpq_mpoly_t A, 
                                           slong i, const fmpq_mpoly_ctx_t ctx)
{
    slong N;

    if ((ulong) i >= (ulong) A->zpoly->length)
    {
        flint_throw(FLINT_ERROR,
                         "Index out of range in fmpz_mpoly_get_term_exp_fmpz");
    }

    N = mpoly_words_per_exp(A->zpoly->bits, ctx->zctx->minfo);
    mpoly_get_monomial_pfmpz(exp, A->zpoly->exps + N*i, A->zpoly->bits,
                                                             ctx->zctx->minfo);
}
