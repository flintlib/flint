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

void gr_mpoly_fit_length_fit_bits(
    gr_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong N = mpoly_words_per_exp(A->bits, mctx);

    if (len > A->coeffs_alloc)
    {
        slong old_alloc = A->coeffs_alloc;
        slong new_alloc = FLINT_MAX(len, 2*old_alloc);
        slong sz = cctx->sizeof_elem;

        A->coeffs_alloc = new_alloc;
        A->coeffs = (gr_ptr) flint_realloc(A->coeffs, new_alloc*sz);
        _gr_vec_init(GR_ENTRY(A->coeffs, old_alloc, sz), new_alloc - old_alloc, cctx);
    }

    if (bits > A->bits)
    {
        slong newN = mpoly_words_per_exp(bits, mctx);
        slong new_exps_alloc = newN*len;
        ulong * t;

        if (len < 1)
        {
            A->bits = bits;
            return;
        }

        t = (ulong *) flint_malloc(new_exps_alloc*sizeof(ulong));

        if (A->length > 0)
            mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length, mctx);

        if (A->exps_alloc > 0)
            flint_free(A->exps);

        A->exps = t;
        A->exps_alloc = new_exps_alloc;
        A->bits = bits;
    }
    else
    {
        if (N*len > A->exps_alloc)
        {
            A->exps_alloc = FLINT_MAX(N*len, 2*A->exps_alloc);
            A->exps = flint_realloc(A->exps, A->exps_alloc*sizeof(ulong));
        }
    }
}
