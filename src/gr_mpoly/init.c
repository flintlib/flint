/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mpoly.h"

void gr_mpoly_init3(
    gr_mpoly_t A,
    slong alloc,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
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
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(MPOLY_MIN_BITS, mctx);
    gr_mpoly_init3(A, alloc, bits, mctx, cctx);
}
