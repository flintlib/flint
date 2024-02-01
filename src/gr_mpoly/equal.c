/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mpoly.h"

/* todo: proper algorithm */
truth_t gr_mpoly_equal(
    const gr_mpoly_t A,
    const gr_mpoly_t B,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    truth_t eq;
    gr_mpoly_t t;

    if (A == B)
        return T_TRUE;

    /* todo: if canonical representation */
    /*
    if (A->length != B->length)
        return T_FALSE;

    ... _gr_vec_equal(A->coeffs, B->coeffs, A->length) ...

    return (0 == mpoly_monomials_cmp(A->exps, A->bits, B->exps, B->bits,
                                                        A->length, ctx->minfo)) ? T_TRUE : T_FALSE;
    */

    gr_mpoly_init(t, mctx, cctx);

    if (gr_mpoly_sub(t, A, B, mctx, cctx) == GR_SUCCESS)
        eq = gr_mpoly_is_zero(t, mctx, cctx);
    else
        eq = T_UNKNOWN;

    gr_mpoly_clear(t, mctx, cctx);

    return eq;
}
