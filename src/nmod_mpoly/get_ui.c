/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

ulong nmod_mpoly_get_ui(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length > WORD(1))
        flint_throw(FLINT_ERROR, "Nonconstant polynomial in nmod_mpoly_get_ui");

    if (A->length == WORD(0))
    {
        return UWORD(0);
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    if (!mpoly_monomial_is_zero(A->exps + N*0, N))
        flint_throw(FLINT_ERROR, "Nonconstant monomial in nmod_mpoly_get_ui");

    return A->coeffs[0];
}
