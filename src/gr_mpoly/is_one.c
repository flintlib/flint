/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

truth_t gr_mpoly_is_one(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)
{
    if (A->length == 0)
        return gr_ctx_is_zero_ring(ctx);

    if (gr_ctx_is_canonical(GR_MPOLY_CCTX(ctx)) == T_TRUE)
    {
        slong N;

        if (A->length != 1)
            return T_FALSE;

        N = mpoly_words_per_exp(A->bits, GR_MPOLY_MCTX(ctx));

        if (!mpoly_monomial_is_zero(A->exps + N*0, N))
            return T_FALSE;

        return gr_is_one(A->coeffs, GR_MPOLY_CCTX(ctx));
    }
    else
    {
        gr_mpoly_t t;
        truth_t res = T_UNKNOWN;
        gr_mpoly_init(t, ctx);
        if (gr_mpoly_one(t, ctx) == GR_SUCCESS)
            res = gr_mpoly_equal(A, t, ctx);
        gr_mpoly_clear(t, ctx);
        return res;
    }
}
