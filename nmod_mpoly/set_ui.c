/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_set_ui(nmod_mpoly_t poly, ulong c, const nmod_mpoly_ctx_t ctx)
{
    mp_limb_t cr;
    slong N;

    NMOD_RED(cr, c, ctx->ffinfo->mod);
    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    if (cr == 0)
    {
        _nmod_mpoly_set_length(poly, 0, ctx);
        return;
    }

    nmod_mpoly_fit_length(poly, 1, ctx);
    poly->coeffs[0] = cr;
    mpoly_monomial_zero(poly->exps, N);
    _nmod_mpoly_set_length(poly, 1, ctx);
}
