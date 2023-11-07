/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mpoly.h"

int gr_mpoly_mul_monomial(
    gr_mpoly_t A,
    const gr_mpoly_t B,
    const gr_mpoly_t C,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    slong i, N, Alen, Blen = B->length;
    ulong ofmask;
    flint_bitcnt_t Abits;
    ulong * Aexps, * Bexps = B->exps, * Cexps = C->exps;
    gr_ptr Ccoeff0;
    int freeCcoeff0 = 0, overflowed = 0;
    int status = GR_SUCCESS;
    TMP_INIT;

    FLINT_ASSERT(C->length == 1);

    if (A == C)
    {
        GR_TMP_INIT(Ccoeff0, cctx);
        status |= gr_set(Ccoeff0, C->coeffs, cctx);
        freeCcoeff0 = 1;
    }
    else
    {
        Ccoeff0 = C->coeffs;
    }

    if (C->exps[0] == 0 && mpoly_monomial_is_zero(C->exps,
                                     mpoly_words_per_exp(C->bits, mctx)))
    {
        status |= gr_mpoly_mul_scalar(A, B, Ccoeff0, mctx, cctx);
        goto cleanup_C;
    }

    TMP_START;

    Abits = FLINT_MAX(B->bits, C->bits);
    N = mpoly_words_per_exp(Abits, mctx);

    if (A == C || Abits != C->bits)
    {
        Cexps = TMP_ARRAY_ALLOC(N, ulong);
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits, 1, mctx);
    }

    if (A == B)
    {
        /* inplace operation on A */
        gr_mpoly_fit_bits(A, Abits, mctx, cctx);
        Bexps = Aexps = A->exps;
    }
    else
    {
        if (Abits != B->bits)
        {
            Bexps = TMP_ARRAY_ALLOC(N*Blen, ulong);
            mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits, Blen, mctx);
        }

        gr_mpoly_fit_length_reset_bits(A, Blen, Abits, mctx, cctx);
        Aexps = A->exps;
    }

    if (Abits > FLINT_BITS)
    {
        for (i = 0; i < Blen; i++)
            mpoly_monomial_add_mp(Aexps + N*i, Bexps + N*i, Cexps + N*0, N);

        for (i = 0; !overflowed && i < Blen; i++)
            overflowed = mpoly_monomial_overflows_mp(Aexps + N*i, N, Abits);
    }
    else
    {
        for (i = 0; i < Blen; i++)
            mpoly_monomial_add(Aexps + N*i, Bexps + N*i, Cexps + N*0, N);

        ofmask = mpoly_overflow_mask_sp(Abits);
        for (i = 0; !overflowed && i < Blen; i++)
            overflowed = mpoly_monomial_overflows(Aexps + N*i, N, ofmask);
    }

    TMP_END;

    /* slightly dirty: repack monomials can handle 1-bit overflown fields */
    if (overflowed)
    {
        ulong * newAexps;
        flint_bitcnt_t newAbits = mpoly_fix_bits(Abits + 1, mctx);
        N = mpoly_words_per_exp(newAbits, mctx);
        newAexps = FLINT_ARRAY_ALLOC(N*A->coeffs_alloc, ulong);
        mpoly_repack_monomials(newAexps, newAbits, A->exps, Abits, Blen, mctx);
        flint_free(A->exps);
        A->exps = newAexps;
        A->bits = newAbits;
        A->exps_alloc = N*A->coeffs_alloc;
    }

/* todo: when we can verify (quickly) that C is invertible */
#if 0
    status |= _gr_vec_mul_scalar(A->coeffs, B->coeffs, Blen, Ccoeff0, cctx);
    _gr_mpoly_set_length(A, Blen, mctx, cctx);
#else
    {
        slong sz = cctx->sizeof_elem;

        Alen = 0;
        for (i = 0; i < Blen; i++)
        {
            mpoly_monomial_set(A->exps + N*Alen, A->exps + N*i, N);
            status |= gr_mul(GR_ENTRY(A->coeffs, Alen, sz), GR_ENTRY(B->coeffs, i, sz), Ccoeff0, cctx);
            Alen += (gr_is_zero(GR_ENTRY(A->coeffs, Alen, sz), cctx) != T_TRUE);
        }

        A->length = Alen;
    }
#endif

cleanup_C:

    if (freeCcoeff0)
        GR_TMP_CLEAR(Ccoeff0, cctx);

    return status;
}
