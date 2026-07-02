/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

/*
    Division by a single term B = c * x^beta.  Every term a*x^alpha of A must
    have x^beta | x^alpha and c | a exactly; the quotient term is then
    (a/c) * x^(alpha - beta).  Since alpha - beta preserves the monomial order
    and is injective, the quotient is already canonical (this also covers the
    scalar case beta = 0).  The coefficient division follows the same rule as
    divides_monagan_pearce: multiply by 1/c when c is a unit, else use gr_div.
*/
static int
_gr_mpoly_divides_monomial(gr_mpoly_t Q,
    const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong sz = cctx->sizeof_elem;
    slong i, N, Alen = A->length;
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

    for (i = 0; i < Alen; i++)
    {
        int divides, cstatus;

        if (exp_bits <= FLINT_BITS)
            divides = mpoly_monomial_divides(T->exps + N*i, Aexps + N*i, Bexps, N, mask);
        else
            divides = mpoly_monomial_divides_mp(T->exps + N*i, Aexps + N*i, Bexps, N, exp_bits);

        if (!divides)
        {
            if (gr_is_zero(GR_ENTRY(A->coeffs, i, sz), cctx) == T_FALSE)
                status = GR_DOMAIN;
            else
                status = GR_UNABLE;
            break;
        }

        if (lc_is_one)
            cstatus = gr_set(GR_ENTRY(T->coeffs, i, sz), GR_ENTRY(A->coeffs, i, sz), cctx);
        else if (lc_is_unit)
            cstatus = gr_mul(GR_ENTRY(T->coeffs, i, sz), GR_ENTRY(A->coeffs, i, sz), inv, cctx);
        else
            cstatus = gr_div(GR_ENTRY(T->coeffs, i, sz), GR_ENTRY(A->coeffs, i, sz), Bc, cctx);

        if (cstatus == GR_DOMAIN)
        {
            status = GR_DOMAIN;
            break;
        }
        else if (cstatus != GR_SUCCESS)
        {
            status |= cstatus;
            break;
        }
    }

    GR_TMP_CLEAR(inv, cctx);

    if (status == GR_SUCCESS)
    {
        _gr_mpoly_set_length(T, Alen, ctx);
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

/*
    TODO: dispatch to gr_mpoly_divides_heap_threaded (port of
    fmpz_mpoly_divides_heap_threaded) for large A when the coefficient ring is
    threadsafe; for now the serial Monagan-Pearce algorithm is used.
*/
int
gr_mpoly_divides(gr_mpoly_t Q,
    const gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        truth_t zero = gr_ctx_is_zero_ring(GR_MPOLY_CCTX(ctx));

        if (zero == T_FALSE)
            return GR_DOMAIN;
        else if (zero == T_TRUE)
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }

    if (A->length == 0)
        return gr_mpoly_zero(Q, ctx);

    /* scalar or monomial divisor: cheaper direct handling */
    if (B->length == 1)
        return _gr_mpoly_divides_monomial(Q, A, B, ctx);

    /* This method is used to overload gr_div, so it should be conservative
       when ctx is something weird. */
    if (gr_ctx_is_integral_domain(GR_MPOLY_CCTX(ctx)) != T_TRUE)
        return GR_UNABLE;

    return gr_mpoly_divides_monagan_pearce(Q, A, B, ctx);
}

