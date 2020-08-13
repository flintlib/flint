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

int fmpz_mpoly_equal_ui(const fmpz_mpoly_t A,
                                           ulong c, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if (c == 0)
        return A->length == 0;

    if (A->length != 1)
        return 0;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (!mpoly_monomial_is_zero(A->exps + N*0, N))
        return 0;

    return fmpz_equal_ui(A->coeffs + 0, c);
}
