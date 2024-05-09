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

/*
    assuming that the exponents are valid and sorted,
    put the polynomial in canonical form
    i.e.
        2*x^e + 3*x^e -> 5x^e
        2*x^e - 2*x^e -> 0
*/
int gr_mpoly_combine_like_terms(
    gr_mpoly_t A,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    gr_method_binary_op add = GR_BINARY_OP(cctx, ADD);
    gr_method_swap_op swap = GR_SWAP_OP(cctx, SWAP);
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(cctx, IS_ZERO);
    int status = GR_SUCCESS;
    slong in, out, N;
    slong sz = cctx->sizeof_elem;

    N = mpoly_words_per_exp(A->bits, mctx);

    out = -WORD(1);

    for (in = 0; in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 &&
                     mpoly_monomial_equal(A->exps + N*out, A->exps + N*in, N))
        {
            status |= add(GR_ENTRY(A->coeffs, out, sz), GR_ENTRY(A->coeffs, out, sz), GR_ENTRY(A->coeffs, in, sz), cctx);
        }
        else
        {
            if (out < 0 || is_zero(GR_ENTRY(A->coeffs, out, sz), cctx) != T_TRUE)
                out++;

            if (out != in)
            {
                mpoly_monomial_set(A->exps + N*out, A->exps + N*in, N);
                swap(GR_ENTRY(A->coeffs, out, sz), GR_ENTRY(A->coeffs, in, sz), cctx);
            }
        }
    }

    if (out < 0 || is_zero(GR_ENTRY(A->coeffs, out, sz), cctx) != T_TRUE)
        out++;

    A->length = out;

    return status;
}
