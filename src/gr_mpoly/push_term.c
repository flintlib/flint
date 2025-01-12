/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

void _gr_mpoly_push_exp_ui(
    gr_mpoly_t A,
    const ulong * exp,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    slong N;
    slong old_length = A->length;
    flint_bitcnt_t exp_bits;

    exp_bits = mpoly_exp_bits_required_ui(exp, mctx);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);
    gr_mpoly_fit_length_fit_bits(A, old_length + 1, exp_bits, ctx);

    A->length = old_length + 1;
    N = mpoly_words_per_exp(A->bits, mctx);
    mpoly_set_monomial_ui(A->exps + N*old_length, exp, A->bits, mctx);
}

int gr_mpoly_push_term_scalar_ui(
    gr_mpoly_t A,
    gr_srcptr c,
    const ulong * exp,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    _gr_mpoly_push_exp_ui(A, exp, ctx);
    return gr_set(GR_ENTRY(A->coeffs, A->length - 1, cctx->sizeof_elem), c, cctx);
}

void _gr_mpoly_push_exp_fmpz(
    gr_mpoly_t A,
    const fmpz * exp,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    slong N;
    slong old_length = A->length;
    flint_bitcnt_t exp_bits;

    exp_bits = mpoly_exp_bits_required_ffmpz(exp, mctx);
    exp_bits = mpoly_fix_bits(exp_bits, mctx);
    gr_mpoly_fit_length_fit_bits(A, old_length + 1, exp_bits, ctx);

    A->length = old_length + 1;
    N = mpoly_words_per_exp(A->bits, mctx);
    mpoly_set_monomial_ffmpz(A->exps + N*old_length, exp, A->bits, mctx);
}

int gr_mpoly_push_term_scalar_fmpz(
    gr_mpoly_t A,
    gr_srcptr c,
    const fmpz * exp,
    gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    _gr_mpoly_push_exp_fmpz(A, exp, ctx);
    return gr_set(GR_ENTRY(A->coeffs, A->length - 1, cctx->sizeof_elem), c, cctx);
}
