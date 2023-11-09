/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_reverse(nmod_mpoly_t A,
                              const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong Blen = B->length;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    if (A != B)
    {
        nmod_mpoly_fit_length_reset_bits(A, Blen, B->bits, ctx);
        A->length = Blen;
        for (i = 0; i < Blen; i++)
            A->coeffs[i] = B->coeffs[Blen - i - 1];
    }
    else
    {
        for (i = 0; i < Blen/2; i++)
            FLINT_SWAP(mp_limb_t, A->coeffs[i], A->coeffs[Blen - i - 1]);
    }

    mpoly_reverse(A->exps, B->exps, Blen, N);
}
