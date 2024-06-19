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

int gr_mpoly_gen(
    gr_mpoly_t A,
    slong var,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    flint_bitcnt_t bits;
    int status = GR_SUCCESS;

    if (var >= mctx->nvars || var < 0)
        return GR_DOMAIN;

    bits = mpoly_gen_bits_required(var, mctx);
    bits = mpoly_fix_bits(bits, mctx);
    gr_mpoly_fit_length_reset_bits(A, 1, bits, mctx, cctx);

    if (bits <= FLINT_BITS)
        mpoly_gen_monomial_sp(A->exps, var, bits, mctx);
    else
        mpoly_gen_monomial_offset_mp(A->exps, var, bits, mctx);

    status |= gr_one(A->coeffs, cctx);
    _gr_mpoly_set_length(A, gr_is_zero(A->coeffs, cctx) != T_TRUE, mctx, cctx);

    return status;
}

truth_t gr_mpoly_is_gen(const gr_mpoly_t A, slong var, const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    truth_t res;

    if (var >= mctx->nvars || mctx->nvars == 0)
        return T_FALSE;

    /* Check for any generator */
    if (var < 0)
        var = -1;

    if (A->length == 1)
    {
        res = mpoly_is_gen(A->exps, var, A->bits, mctx) ? T_TRUE : T_FALSE;

        if (res != T_TRUE)
            return res;

        res = gr_is_one(A->coeffs, cctx);
    }
    else
    {
        /* todo: cheaper check when possible */
        gr_mpoly_t t;

        gr_mpoly_init(t, mctx, cctx);

        if (gr_mpoly_gen(t, var, mctx, cctx) != GR_SUCCESS)
            res = T_UNKNOWN;
        else
            res = gr_mpoly_equal(A, t, mctx, cctx);

        gr_mpoly_clear(t, mctx, cctx);
    }

    return res;
}
