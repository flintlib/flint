/*
    Copyright (C) 2020, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void _fmpz_mod_mpoly_push_exp_ffmpz(
    fmpz_mod_mpoly_t A,
    const fmpz * exp,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;
    slong old_length = A->length;
    flint_bitcnt_t exp_bits;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mod_mpoly_fit_length_fit_bits(A, old_length + 1, exp_bits, ctx);

    A->length = old_length + 1;
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_set_monomial_ffmpz(A->exps + N*old_length, exp, A->bits, ctx->minfo);
}

void _fmpz_mod_mpoly_push_exp_pfmpz(
    fmpz_mod_mpoly_t A,
    fmpz * const * exp,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;
    slong old_length = A->length;
    flint_bitcnt_t exp_bits;

    exp_bits = mpoly_exp_bits_required_pfmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mod_mpoly_fit_length_fit_bits(A, old_length + 1, exp_bits, ctx);

    A->length = old_length + 1;
    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    mpoly_set_monomial_pfmpz(A->exps + N*old_length, exp, A->bits, ctx->minfo);
}


void fmpz_mod_mpoly_push_term_fmpz_fmpz(
    fmpz_mod_mpoly_t A,
    const fmpz_t c,
    fmpz * const * exp,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_push_exp_pfmpz(A, exp, ctx);
    fmpz_mod_set_fmpz(A->coeffs + A->length - 1, c, ctx->ffinfo);
}

void fmpz_mod_mpoly_push_term_ui_fmpz(
    fmpz_mod_mpoly_t A,
    ulong c,
    fmpz * const * exp,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_push_exp_pfmpz(A, exp, ctx);
    fmpz_mod_set_ui(A->coeffs + A->length - 1, c, ctx->ffinfo);
}

void fmpz_mod_mpoly_push_term_si_fmpz(
    fmpz_mod_mpoly_t A,
    slong c,
    fmpz * const * exp,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_push_exp_pfmpz(A, exp, ctx);
    fmpz_mod_set_si(A->coeffs + A->length - 1, c, ctx->ffinfo);
}

