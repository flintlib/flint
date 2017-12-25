/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int nmod_mpoly_equal_ui(const nmod_mpoly_t poly,
                                           ulong c, const nmod_mpoly_ctx_t ctx)
{
    slong N;
    mp_limb_t cr;

    NMOD_RED(cr, c, ctx->ffinfo->mod);

    if (cr == 0)
        return poly->length == 0;

    if (poly->length != 1)
        return 0;

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    if (!mpoly_monomial_is_zero(poly->exps + N*0, N))
        return 0;

    return poly->coeffs[0] == cr;
}
