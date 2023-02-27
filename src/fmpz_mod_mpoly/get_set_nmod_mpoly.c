/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


void _fmpz_mod_mpoly_set_nmod_mpoly(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx,
    const nmod_mpoly_t nA,
    const nmod_mpoly_ctx_t nctx)
{
    slong i, N = mpoly_words_per_exp(nA->bits, ctx->minfo);
    fmpz_mod_mpoly_fit_length_reset_bits(A, nA->length, nA->bits, ctx);
    mpoly_copy_monomials(A->exps, nA->exps, nA->length, N);
    for (i = 0; i < nA->length; i++)
        fmpz_set_ui(A->coeffs + i, nA->coeffs[i]);
    _fmpz_mod_mpoly_set_length(A, nA->length, ctx);
}


void _fmpz_mod_mpoly_get_nmod_mpoly(
    nmod_mpoly_t nA,
    const nmod_mpoly_ctx_t nctx,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N = mpoly_words_per_exp(A->bits, ctx->minfo);
    nmod_mpoly_fit_length_reset_bits(nA, A->length, A->bits, nctx);
    mpoly_copy_monomials(nA->exps, A->exps, A->length, N);
    for (i = 0; i < A->length; i++)
        nA->coeffs[i] = fmpz_get_ui(A->coeffs + i);
    nA->length = A->length;
}
