/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "gr_mpoly.h"

/*
    Square of a polynomial with at most one term: (a X^e)^2 = a^2 X^{2e}.
    This is a specialization of gr_mpoly_mul_monomial for the case B == C.
    Only the exponent doubling and a single coefficient square are needed, and
    the operation is performed in place when A aliases B.
*/

int gr_mpoly_sqr_monomial(
    gr_mpoly_t A,
    const gr_mpoly_t B,
    gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong N;
    flint_bitcnt_t Abits;
    ulong ofmask;
    ulong * Aexps, * Bexps = B->exps;
    int overflowed = 0;
    int status = GR_SUCCESS;

    FLINT_ASSERT(B->length <= 1);

    if (B->length == 0)
        return gr_mpoly_zero(A, ctx);

    Abits = B->bits;
    N = mpoly_words_per_exp(Abits, mctx);

    if (A == B)
    {
        /* in place: exps and coeff are updated where they sit */
        Aexps = A->exps;
    }
    else
    {
        gr_mpoly_fit_length_reset_bits(A, 1, Abits, ctx);
        Aexps = A->exps;
    }

    /* double the exponent vector: 2e = e + e */
    if (Abits > FLINT_BITS)
    {
        mpoly_monomial_add_mp(Aexps, Bexps, Bexps, N);
        overflowed = mpoly_monomial_overflows_mp(Aexps, N, Abits);
    }
    else
    {
        mpoly_monomial_add(Aexps, Bexps, Bexps, N);
        ofmask = mpoly_overflow_mask_sp(Abits);
        overflowed = mpoly_monomial_overflows(Aexps, N, ofmask);
    }

    /* slightly dirty: repack monomials can handle 1-bit overflown fields */
    if (overflowed)
    {
        ulong * newAexps;
        flint_bitcnt_t newAbits = mpoly_fix_bits(Abits + 1, mctx);
        N = mpoly_words_per_exp(newAbits, mctx);
        newAexps = FLINT_ARRAY_ALLOC(N*A->coeffs_alloc, ulong);
        mpoly_repack_monomials(newAexps, newAbits, A->exps, Abits, 1, mctx);
        flint_free(A->exps);
        A->exps = newAexps;
        A->bits = newAbits;
        A->exps_alloc = N*A->coeffs_alloc;
    }

    /* square the coefficient (gr_sqr is aliasing-safe, so the A == B case is
       handled directly) */
    status |= gr_sqr(A->coeffs, B->coeffs, cctx);
    A->length = (gr_is_zero(A->coeffs, cctx) != T_TRUE);

    return status;
}
