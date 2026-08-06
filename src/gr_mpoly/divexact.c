/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2017 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"
#include "threaded_divmod.h" /* for the exposed _gr_mpoly_divrem_mp kernel */

static int
_gr_mpoly_divexact_monomial(gr_mpoly_t Q,
    const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    slong i, qlen, N, Alen = A->length;
    flint_bitcnt_t exp_bits;
    ulong * Aexps = A->exps, * Bexps = B->exps;
    gr_srcptr Bc = B->coeffs;
    int freeA = 0, freeB = 0, lc_is_unit, lc_is_one;
    ulong mask;
    gr_ptr inv;
    int status = GR_SUCCESS;
    gr_mpoly_t T;

    FLINT_ASSERT(B->length == 1);

    if (gr_is_zero(B->coeffs, cctx) != T_FALSE)
        return GR_UNABLE;

    exp_bits = FLINT_MAX(A->bits, B->bits);
    N = mpoly_words_per_exp(exp_bits, mctx);

    if (exp_bits > A->bits)
    {
        Aexps = (ulong *) flint_malloc(N*Alen*sizeof(ulong));
        mpoly_repack_monomials(Aexps, exp_bits, A->exps, A->bits, Alen, mctx);
        freeA = 1;
    }

    if (exp_bits > B->bits)
    {
        Bexps = (ulong *) flint_malloc(N*sizeof(ulong));
        mpoly_repack_monomials(Bexps, exp_bits, B->exps, B->bits, 1, mctx);
        freeB = 1;
    }

    /* mask with the high bit of each exponent field set (single-word fields) */
    mask = 0;
    if (exp_bits <= FLINT_BITS)
    {
        flint_bitcnt_t j;
        for (j = 0; j < FLINT_BITS/exp_bits; j++)
            mask = (mask << exp_bits) + (UWORD(1) << (exp_bits - 1));
    }

    gr_mpoly_init3(T, Alen, exp_bits, ctx);

    GR_TMP_INIT(inv, cctx);
    lc_is_one = (gr_is_one(Bc, cctx) == T_TRUE);
    lc_is_unit = lc_is_one || (gr_inv(inv, Bc, cctx) == GR_SUCCESS);

    qlen = 0;
    for (i = 0; i < Alen; i++)
    {
        int divides, cstatus;

        divides = mpoly_monomial_divides_any_bits(T->exps + N*qlen, Aexps + N*i, Bexps, N, mask, exp_bits);

        /* Since the division is known to be exact, spurious terms must correspond
           to inexact/unproved representation of actual zeros, which are thus
           safe to discard (e.g. we must have [+/-0.1]*x / y -> 0). */
        if (!divides)
            continue;

        if (lc_is_one)
            cstatus = gr_set(GR_ENTRY(T->coeffs, qlen, sz), GR_ENTRY(A->coeffs, i, sz), cctx);
        else if (lc_is_unit)
            cstatus = gr_mul(GR_ENTRY(T->coeffs, qlen, sz), GR_ENTRY(A->coeffs, i, sz), inv, cctx);
        else
            cstatus = gr_divexact(GR_ENTRY(T->coeffs, qlen, sz), GR_ENTRY(A->coeffs, i, sz), Bc, cctx);

        if (cstatus != GR_SUCCESS)
        {
            status |= cstatus;
            break;
        }

        qlen++;
    }

    GR_TMP_CLEAR(inv, cctx);

    if (status == GR_SUCCESS)
    {
        _gr_mpoly_set_length(T, qlen, ctx);
        gr_mpoly_swap(Q, T, ctx);
    }
    else
    {
        GR_IGNORE(gr_mpoly_zero(Q, ctx));
    }

    gr_mpoly_clear(T, ctx);

    if (freeA)
        flint_free(Aexps);
    if (freeB)
        flint_free(Bexps);

    return status;
}


int
gr_mpoly_divexact(gr_mpoly_t Q,
    const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);

    if (B->length == 0)
    {
        truth_t zero = gr_ctx_is_zero_ring(cctx);

        if (zero == T_FALSE)
            return GR_DOMAIN;
        else if (zero == T_TRUE)
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }

    if (A->length == 0)
        return gr_mpoly_zero(Q, ctx);

    if (B->length == 1)
        return _gr_mpoly_divexact_monomial(Q, A, B, ctx);

    if (A->length > 500 && B->length > 2 &&
            gr_ctx_is_threadsafe(cctx) == T_TRUE &&
            flint_get_num_available_threads() > 1)
        return gr_mpoly_divexact_heap_threaded(Q, A, B, ctx);

    return _gr_mpoly_divrem_mp(Q, NULL, A, B, 0, 1, ctx);
}
