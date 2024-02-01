/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mpoly.h"

int gr_mpoly_randtest_bound(gr_mpoly_t A, flint_rand_t state,
                 slong length, ulong exp_bound, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong i, j, nvars = mctx->nvars;
    ulong * exp;
    int status = GR_SUCCESS;
    TMP_INIT;

    TMP_START;
    exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    status = gr_mpoly_zero(A, mctx, cctx);

    gr_mpoly_fit_length_reset_bits(A, 0, MPOLY_MIN_BITS, mctx, cctx);

    for (i = 0; i < length; i++)
    {
        for (j = 0; j < nvars; j++)
            exp[j] = n_randint(state, exp_bound);
        _gr_mpoly_push_exp_ui(A, exp, mctx, cctx);
        status |= gr_randtest(GR_ENTRY(A->coeffs, A->length - 1, cctx->sizeof_elem), state, cctx);
    }

    gr_mpoly_sort_terms(A, mctx, cctx);
    status |= gr_mpoly_combine_like_terms(A, mctx, cctx);

    TMP_END;

    return status;
}
