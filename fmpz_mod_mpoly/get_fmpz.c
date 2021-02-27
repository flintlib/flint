/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_get_fmpz(
    fmpz_t c,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length > 1)
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_get_fmpz: nonconstant polynomial");

    if (A->length < 1)
    {
        fmpz_zero(c);
        return;
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    if (!mpoly_monomial_is_zero(A->exps + N*0, N))
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_get_fmpz: nonconstant polynomial");

    fmpz_set(c, A->coeffs + 0);
}
