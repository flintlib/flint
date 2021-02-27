/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_set(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    if (A == B)
        return;

    fmpz_mod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);

    _fmpz_vec_set(A->coeffs, B->coeffs, B->length);
    mpoly_copy_monomials(A->exps, B->exps, B->length, N);
    _fmpz_mod_mpoly_set_length(A, B->length, ctx);
}
