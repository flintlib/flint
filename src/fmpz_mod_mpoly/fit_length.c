/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


void fmpz_mod_mpoly_fit_length(
    fmpz_mod_mpoly_t A,
    slong len,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    _fmpz_mod_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc,
                               &A->exps, &A->exps_alloc, N, len);
}
