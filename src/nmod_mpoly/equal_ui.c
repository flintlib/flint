/*
    Copyright (C) 2017, 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int nmod_mpoly_equal_ui(const nmod_mpoly_t A,
                                           ulong c, const nmod_mpoly_ctx_t ctx)
{
    slong N;

    if (c >= ctx->mod.n)
        NMOD_RED(c, c, ctx->mod);

    if (c == 0)
        return A->length == 0;

    if (A->length != 1)
        return 0;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (!mpoly_monomial_is_zero(A->exps + N*0, N))
        return 0;

    return A->coeffs[0] == c;
}
