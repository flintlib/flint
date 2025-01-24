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

void gr_mpoly_init3(
    gr_mpoly_t A,
    slong alloc,
    flint_bitcnt_t bits,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong N = mpoly_words_per_exp(bits, mctx);

    if (alloc > 0)
    {
        A->coeffs_alloc = alloc;
        A->coeffs = (gr_ptr) flint_malloc(alloc * cctx->sizeof_elem);
        _gr_vec_init(A->coeffs, alloc, cctx);
        A->exps_alloc = N * alloc;
        A->exps = FLINT_ARRAY_ALLOC(A->exps_alloc, ulong);
    }
    else
    {
        A->coeffs = NULL;
        A->exps = NULL;
        A->coeffs_alloc = 0;
        A->exps_alloc = 0;
    }

    A->length = 0;
    A->bits = bits;
}

void gr_mpoly_init2(
    gr_mpoly_t A,
    slong alloc,
    gr_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(MPOLY_MIN_BITS, GR_MPOLY_MCTX(ctx));
    gr_mpoly_init3(A, alloc, bits, ctx);
}
