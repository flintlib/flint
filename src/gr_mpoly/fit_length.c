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

void _gr_mpoly_fit_length(
    gr_ptr * coeffs,
    slong * coeffs_alloc,
    ulong ** exps,
    slong * exps_alloc,
    slong N,
    slong length,
    gr_mpoly_ctx_t ctx)
{
    if (length > *coeffs_alloc)
    {
        gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);

        slong old_alloc = *coeffs_alloc;
        slong new_alloc = FLINT_MAX(length, old_alloc*2);
        slong sz = cctx->sizeof_elem;

        *coeffs_alloc = new_alloc;
        *coeffs = (gr_ptr) flint_realloc(*coeffs, new_alloc*sz);
        _gr_vec_init(GR_ENTRY(*coeffs, old_alloc, sz), new_alloc - old_alloc, cctx);
    }

    if (N*length > *exps_alloc)
    {
        *exps_alloc = FLINT_MAX(N*length, *exps_alloc*2);
        *exps = (ulong *) flint_realloc(*exps, *exps_alloc*sizeof(ulong));
    }
}

void gr_mpoly_fit_length(
    gr_mpoly_t A,
    slong len,
    gr_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, GR_MPOLY_MCTX(ctx));

    _gr_mpoly_fit_length(&A->coeffs, &A->coeffs_alloc,
                               &A->exps, &A->exps_alloc, N, len, ctx);
}
