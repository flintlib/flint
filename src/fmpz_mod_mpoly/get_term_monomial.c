/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_get_term_monomial(
    fmpz_mod_mpoly_t M,
    const fmpz_mod_mpoly_t A,
    slong i,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    if (i >= (ulong) A->length)
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_get_term_monomial: index out of range");

    fmpz_mod_mpoly_fit_length_reset_bits(M, 1, bits, ctx);

    mpoly_monomial_set(M->exps + N*0, A->exps + N*i, N);
    fmpz_one(M->coeffs + 0);
    _fmpz_mod_mpoly_set_length(M, 1, ctx);
}
