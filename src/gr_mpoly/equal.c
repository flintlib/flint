/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

truth_t gr_mpoly_equal(
    const gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    int canonical;

    if (A == B)
        return T_TRUE;

    canonical = 1;

    if (gr_ctx_is_canonical(cctx) != T_TRUE)
    {
        slong i, sz = cctx->sizeof_elem;

        for (i = 0; canonical && i < A->length; i++)
            if (gr_is_zero(GR_ENTRY(A->coeffs, i, sz), cctx) != T_FALSE)
                canonical = 0;

        for (i = 0; canonical && i < B->length; i++)
            if (gr_is_zero(GR_ENTRY(B->coeffs, i, sz), cctx) != T_FALSE)
                canonical = 0;
    }

    if (canonical)
    {
        if (A->length != B->length)
            return T_FALSE;

        if (0 != mpoly_monomials_cmp(A->exps, A->bits, B->exps, B->bits, A->length, mctx))
            return T_FALSE;

        return _gr_vec_equal(A->coeffs, B->coeffs, A->length, cctx);
    }
    else
    {
        /* todo: a better fallback algorithm */
        truth_t eq;
        gr_mpoly_t t;
        gr_mpoly_init(t, ctx);

        if (gr_mpoly_sub(t, A, B, ctx) == GR_SUCCESS)
            eq = gr_mpoly_is_zero(t, ctx);
        else
            eq = T_UNKNOWN;

        gr_mpoly_clear(t, ctx);

        return eq;
    }
}
