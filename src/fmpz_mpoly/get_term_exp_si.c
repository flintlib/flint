/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_get_term_exp_si(slong * exp, const fmpz_mpoly_t A,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR,
                           "Index out of range in fmpz_mpoly_get_term_exp_si");
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_get_monomial_si(exp, A->exps + N*i, A->bits, ctx->minfo);
}
