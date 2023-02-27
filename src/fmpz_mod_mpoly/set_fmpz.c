/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_set_fmpz_mod(
    fmpz_mod_mpoly_t A,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx->ffinfo));

    fmpz_mod_mpoly_fit_length(A, 1, ctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, ctx->minfo));
    fmpz_set(A->coeffs + 0, c);
    _fmpz_mod_mpoly_set_length(A, !fmpz_is_zero(A->coeffs + 0), ctx);
}


void fmpz_mod_mpoly_set_fmpz(
    fmpz_mod_mpoly_t A,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_fit_length(A, 1, ctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, ctx->minfo));
    fmpz_mod_set_fmpz(A->coeffs + 0, c, ctx->ffinfo);
    _fmpz_mod_mpoly_set_length(A, !fmpz_is_zero(A->coeffs + 0), ctx);
}


void fmpz_mod_mpoly_set_ui(
    fmpz_mod_mpoly_t A,
    slong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_fit_length(A, 1, ctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, ctx->minfo));
    fmpz_mod_set_ui(A->coeffs + 0, c, ctx->ffinfo);
    _fmpz_mod_mpoly_set_length(A, !fmpz_is_zero(A->coeffs + 0), ctx);
}


void fmpz_mod_mpoly_set_si(
    fmpz_mod_mpoly_t A,
    slong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_fit_length(A, 1, ctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, ctx->minfo));
    fmpz_mod_set_si(A->coeffs + 0, c, ctx->ffinfo);
    _fmpz_mod_mpoly_set_length(A, !fmpz_is_zero(A->coeffs + 0), ctx);
}
