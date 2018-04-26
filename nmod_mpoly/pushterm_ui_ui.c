/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/*
    emplaceterm assumes that c is valid modulo ctx->ffinfo->mod.n
*/

void _nmod_mpoly_emplacebackterm_ui_ui(nmod_mpoly_t poly,
                    mp_limb_t c, const ulong * exp, const nmod_mpoly_ctx_t ctx)
{
    mp_bitcnt_t exp_bits;
    slong N;
    FLINT_ASSERT(c < ctx->ffinfo->mod.n);

    exp_bits = mpoly_exp_bits_required_ui(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    nmod_mpoly_fit_bits(poly, exp_bits, ctx);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);

    nmod_mpoly_fit_length(poly, poly->length + 1, ctx);
    poly->coeffs[poly->length] = c;
    mpoly_set_monomial_ui(poly->exps + N*poly->length, exp, poly->bits, ctx->minfo);
    poly->length++; /* safe because length is increasing */
}


void nmod_mpoly_pushterm_ui_ui(nmod_mpoly_t poly,
                        ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx)
{
    ulong C;
    NMOD_RED(C, c, ctx->ffinfo->mod);
    _nmod_mpoly_emplacebackterm_ui_ui(poly, C, exp, ctx);    
}
