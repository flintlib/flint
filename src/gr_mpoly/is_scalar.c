/*
    Copyright (C) 2026 Oscar Benjamin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

truth_t
gr_mpoly_is_scalar(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    truth_t res = T_TRUE;
    slong i, N;

    if (gr_ctx_is_canonical(cctx) == T_TRUE)
    {
        N = mpoly_words_per_exp(A->bits, mctx);

        for (i = 0; i < A->length; i++)
            if (!mpoly_monomial_is_zero(A->exps + N*i, N))
                return T_FALSE;

        return T_TRUE;
    }

    N = mpoly_words_per_exp(A->bits, mctx);

    for (i = 0; i < A->length; i++)
    {
        if (!mpoly_monomial_is_zero(A->exps + N*i, N))
        {
            truth_t is_zero = gr_is_zero(GR_ENTRY(A->coeffs, i, cctx->sizeof_elem), cctx);

            if (is_zero == T_FALSE)
                return T_FALSE;
            if (is_zero == T_UNKNOWN)
                res = T_UNKNOWN;
        }
    }

    return res;
}
