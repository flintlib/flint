/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_gen(nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;

    bits = mpoly_gen_bits_required(var, ctx->minfo);
    bits = mpoly_fix_bits(bits, ctx->minfo);

    nmod_mpoly_fit_length(A, WORD(1), ctx);
    nmod_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    A->coeffs[0] = UWORD(1);
    if (bits <= FLINT_BITS)
        mpoly_gen_monomial_sp(A->exps, var, bits, ctx->minfo);
    else
        mpoly_gen_monomial_offset_mp(A->exps, var, bits, ctx->minfo);

    _nmod_mpoly_set_length(A, WORD(1), ctx);
}
