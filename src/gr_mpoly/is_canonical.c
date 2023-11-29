/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mpoly.h"

truth_t gr_mpoly_is_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong N = mpoly_words_per_exp(A->bits, mctx);
    slong i;
    truth_t ok;

    if (A->length > A->coeffs_alloc)
        return T_FALSE;

    if (N * A->length > A->exps_alloc)
        return T_FALSE;

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, mctx))
        return T_FALSE;

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, mctx))
        return T_FALSE;

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, mctx))
        return T_FALSE;

    ok = T_TRUE;

    for (i = 0; i < A->length; i++)
        ok = truth_and(ok, truth_not(gr_is_zero(GR_ENTRY(A->coeffs, i, cctx->sizeof_elem), cctx)));

    return ok;
}

void gr_mpoly_assert_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong N = mpoly_words_per_exp(A->bits, mctx);
    slong i;

    if (A->length > A->coeffs_alloc)
        flint_throw(FLINT_ERROR, "Polynomial coefficient allocation is bad");

    if (N * A->length > A->exps_alloc)
        flint_throw(FLINT_ERROR, "Polynomial exponent allocation is bad");

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, mctx))
        flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, mctx))
        flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, mctx))
        flint_throw(FLINT_ERROR, "Polynomial exponents out of order");

    for (i = 0; i < A->length; i++)
    {
        if (gr_is_zero(GR_ENTRY(A->coeffs, i, cctx->sizeof_elem), cctx) == T_TRUE)
            flint_throw(FLINT_ERROR, "Polynomial has a zero coefficient");
    }
}
