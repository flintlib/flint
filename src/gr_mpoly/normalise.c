/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

/*
    remove zero coefficients
*/
void
_gr_mpoly_normalise(gr_mpoly_t A, gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    gr_method_swap_op swap = GR_SWAP_OP(cctx, SWAP);
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(cctx, IS_ZERO);
    slong in, out, N;
    slong sz = cctx->sizeof_elem;

    N = mpoly_words_per_exp(A->bits, mctx);

    out = -WORD(1);

    for (in = 0; in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out < 0 || is_zero(GR_ENTRY(A->coeffs, out, sz), cctx) != T_TRUE)
            out++;

        if (out != in)
        {
            mpoly_monomial_set(A->exps + N*out, A->exps + N*in, N);
            swap(GR_ENTRY(A->coeffs, out, sz), GR_ENTRY(A->coeffs, in, sz), cctx);
        }
    }

    if (out < 0 || is_zero(GR_ENTRY(A->coeffs, out, sz), cctx) != T_TRUE)
        out++;

    A->length = out;
}
