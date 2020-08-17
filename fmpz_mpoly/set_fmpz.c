/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_set_fmpz(fmpz_mpoly_t A,
                                    const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (fmpz_is_zero(c))
    {
        _fmpz_mpoly_set_length(A, 0, ctx);
        return;
    }

    fmpz_mpoly_fit_length(A, 1, ctx);
    fmpz_set(A->coeffs + 0, c);
    mpoly_monomial_zero(A->exps + N*0, N);
    _fmpz_mpoly_set_length(A, 1, ctx);
}
