/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "gr_mpoly.h"

int gr_mpoly_set_scalar(gr_mpoly_t A,
    gr_srcptr c,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;

    if (gr_is_zero(c, cctx) == T_TRUE)
        return gr_mpoly_zero(A, mctx, cctx);

    gr_mpoly_fit_length(A, 1, mctx, cctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set(A->coeffs, c, cctx);
    _gr_mpoly_set_length(A, 1, mctx, cctx);

    return status;
}

int gr_mpoly_set_ui(gr_mpoly_t A,
    ulong c,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;

    if (c == 0)
        return gr_mpoly_zero(A, mctx, cctx);

    gr_mpoly_fit_length(A, 1, mctx, cctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set_ui(A->coeffs, c, cctx);

    if (gr_is_zero(A->coeffs, cctx) == T_TRUE)
        _gr_mpoly_set_length(A, 0, mctx, cctx);
    else
        _gr_mpoly_set_length(A, 1, mctx, cctx);

    return status;
}

int gr_mpoly_set_si(gr_mpoly_t A,
    slong c,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;

    if (c == 0)
        return gr_mpoly_zero(A, mctx, cctx);

    gr_mpoly_fit_length(A, 1, mctx, cctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set_si(A->coeffs, c, cctx);

    if (gr_is_zero(A->coeffs, cctx) == T_TRUE)
        _gr_mpoly_set_length(A, 0, mctx, cctx);
    else
        _gr_mpoly_set_length(A, 1, mctx, cctx);

    return status;
}

int gr_mpoly_set_fmpz(gr_mpoly_t A,
    const fmpz_t c,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;

    if (fmpz_is_zero(c))
        return gr_mpoly_zero(A, mctx, cctx);

    gr_mpoly_fit_length(A, 1, mctx, cctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set_fmpz(A->coeffs, c, cctx);

    if (gr_is_zero(A->coeffs, cctx) == T_TRUE)
        _gr_mpoly_set_length(A, 0, mctx, cctx);
    else
        _gr_mpoly_set_length(A, 1, mctx, cctx);

    return status;
}

int gr_mpoly_set_fmpq(gr_mpoly_t A,
    const fmpq_t c,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    int status;

    if (fmpq_is_zero(c))
        return gr_mpoly_zero(A, mctx, cctx);

    gr_mpoly_fit_length(A, 1, mctx, cctx);
    mpoly_monomial_zero(A->exps, mpoly_words_per_exp(A->bits, mctx));

    status = gr_set_fmpq(A->coeffs, c, cctx);

    if (gr_is_zero(A->coeffs, cctx) == T_TRUE)
        _gr_mpoly_set_length(A, 0, mctx, cctx);
    else
        _gr_mpoly_set_length(A, 1, mctx, cctx);

    return status;
}
