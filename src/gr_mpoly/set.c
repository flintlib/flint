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

int gr_mpoly_set(
    gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong N;
    int status = GR_SUCCESS;

    if (A == B)
        return GR_SUCCESS;

    N = mpoly_words_per_exp(B->bits, mctx);

    gr_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);

    status = _gr_vec_set(A->coeffs, B->coeffs, B->length, cctx);
    mpoly_copy_monomials(A->exps, B->exps, B->length, N);
    _gr_mpoly_set_length(A, B->length, ctx);

    return status;
}
