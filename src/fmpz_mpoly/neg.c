/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_neg(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if (A != B)
    {
        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        fmpz_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
        mpn_copyi(A->exps, B->exps, N*B->length);
    }
    _fmpz_vec_neg(A->coeffs, B->coeffs, B->length);
    _fmpz_mpoly_set_length(A, B->length, ctx);
}
