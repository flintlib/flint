/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_neg(nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    if (A != B)
    {
        slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
        nmod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
        mpoly_copy_monomials(A->exps, B->exps, B->length, N);
        _nmod_mpoly_set_length(A, B->length, ctx);
    }

    _nmod_vec_neg(A->coeffs, B->coeffs, B->length, ctx->mod);
}
