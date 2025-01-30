/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "mpoly.h"
#include "gr_mpoly.h"

int gr_mpoly_set_scalar(gr_mpoly_t A,
    gr_srcptr c,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    int status;

    if (gr_is_zero(c, cctx) == T_TRUE)
        return gr_mpoly_zero(A, ctx);

    gr_mpoly_fit_length(A, 1, ctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set(A->coeffs, c, cctx);
    _gr_mpoly_set_length(A, 1, ctx);

    return status;
}

int gr_mpoly_set_ui(gr_mpoly_t A,
    ulong c,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    int status;

    if (c == 0)
        return gr_mpoly_zero(A, ctx);

    gr_mpoly_fit_length(A, 1, ctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set_ui(A->coeffs, c, cctx);

    if (gr_is_zero(A->coeffs, cctx) == T_TRUE)
        _gr_mpoly_set_length(A, 0, ctx);
    else
        _gr_mpoly_set_length(A, 1, ctx);

    return status;
}

int gr_mpoly_set_si(gr_mpoly_t A,
    slong c,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    int status;

    if (c == 0)
        return gr_mpoly_zero(A, ctx);

    gr_mpoly_fit_length(A, 1, ctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set_si(A->coeffs, c, cctx);

    if (gr_is_zero(A->coeffs, cctx) == T_TRUE)
        _gr_mpoly_set_length(A, 0, ctx);
    else
        _gr_mpoly_set_length(A, 1, ctx);

    return status;
}

int gr_mpoly_set_fmpz(gr_mpoly_t A,
    const fmpz_t c,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    int status;

    if (fmpz_is_zero(c))
        return gr_mpoly_zero(A, ctx);

    gr_mpoly_fit_length(A, 1, ctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set_fmpz(A->coeffs, c, cctx);

    if (gr_is_zero(A->coeffs, cctx) == T_TRUE)
        _gr_mpoly_set_length(A, 0, ctx);
    else
        _gr_mpoly_set_length(A, 1, ctx);

    return status;
}

int gr_mpoly_set_fmpq(gr_mpoly_t A,
    const fmpq_t c,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    int status;

    if (fmpq_is_zero(c))
        return gr_mpoly_zero(A, ctx);

    gr_mpoly_fit_length(A, 1, ctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set_fmpq(A->coeffs, c, cctx);

    if (gr_is_zero(A->coeffs, cctx) == T_TRUE)
        _gr_mpoly_set_length(A, 0, ctx);
    else
        _gr_mpoly_set_length(A, 1, ctx);

    return status;
}
