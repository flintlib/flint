/*
    Copyright (C) 2018, 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void _fmpz_mpoly_push_exp_ffmpz(fmpz_mpoly_t A,
                                 const fmpz * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong N;
    slong old_length = A->length;
    flint_bitcnt_t exp_bits;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(A, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fmpz_mpoly_fit_length(A, old_length + 1, ctx);
    A->length = old_length + 1;
    mpoly_set_monomial_ffmpz(A->exps + N*old_length, exp, A->bits, ctx->minfo);
}

void _fmpz_mpoly_push_exp_pfmpz(fmpz_mpoly_t A,
                                fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong N;
    slong old_length = A->length;
    flint_bitcnt_t exp_bits;

    exp_bits = mpoly_exp_bits_required_pfmpz(exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(A, exp_bits, ctx);

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fmpz_mpoly_fit_length(A, old_length + 1, ctx);
    A->length = old_length + 1;
    mpoly_set_monomial_pfmpz(A->exps + N*old_length, exp, A->bits, ctx->minfo);
}


void fmpz_mpoly_push_term_fmpz_fmpz(fmpz_mpoly_t A,
                const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_push_exp_pfmpz(A, exp, ctx);
    fmpz_set(A->coeffs + A->length - 1, c);
}

void fmpz_mpoly_push_term_ui_fmpz(fmpz_mpoly_t A,
                       ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_push_exp_pfmpz(A, exp, ctx);
    fmpz_set_ui(A->coeffs + A->length - 1, c);
}

void fmpz_mpoly_push_term_si_fmpz(fmpz_mpoly_t A,
                       slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_push_exp_pfmpz(A, exp, ctx);
    fmpz_set_si(A->coeffs + A->length - 1, c);
}
