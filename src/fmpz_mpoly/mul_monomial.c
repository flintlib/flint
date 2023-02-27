/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/* TODO: decouple the exp/coeff alloc in fmpz_mpoly and move this to mpoly */
void fmpz_mpoly_mul_monomial(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, Blen = B->length;
    ulong ofmask;
    flint_bitcnt_t Abits;
    ulong * Aexps, * Bexps = B->exps, * Cexps = C->exps;
    fmpz Ccoeff0 = C->coeffs[0];
    int freeCcoeff0 = 0, overflowed = 0;
    TMP_INIT;

    FLINT_ASSERT(C->length == 1);

    if (A == C)
    {
        freeCcoeff0 = 1;
        fmpz_init_set(&Ccoeff0, C->coeffs + 0);
    }

    if (C->exps[0] == 0 && mpoly_monomial_is_zero(C->exps,
                                     mpoly_words_per_exp(C->bits, ctx->minfo)))
    {
        fmpz_mpoly_scalar_mul_fmpz(A, B, &Ccoeff0, ctx);
        goto cleanup_C;
    }

    TMP_START;

    Abits = FLINT_MAX(B->bits, C->bits);
    N = mpoly_words_per_exp(Abits, ctx->minfo);

    if (A == C || Abits != C->bits)
    {
        Cexps = TMP_ARRAY_ALLOC(N, ulong);
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits, 1, ctx->minfo);
    }

    if (A == B)
    {
        /* inplace operation on A */
        fmpz_mpoly_fit_bits(A, Abits, ctx);
        Bexps = Aexps = A->exps;
    }
    else
    {
        if (Abits != B->bits)
        {
            Bexps = TMP_ARRAY_ALLOC(N*Blen, ulong);
            mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits, Blen, ctx->minfo);
        }

        fmpz_mpoly_fit_length_reset_bits(A, Blen, Abits, ctx);
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

    /* slightly dirty: repack monomials can handle 1-bit overfown fields */
    if (overflowed)
    {
        ulong * newAexps;
        flint_bitcnt_t newAbits = mpoly_fix_bits(Abits + 1, ctx->minfo);
        N = mpoly_words_per_exp(newAbits, ctx->minfo);
        newAexps = FLINT_ARRAY_ALLOC(N*A->alloc, ulong);
        mpoly_repack_monomials(newAexps, newAbits, A->exps, Abits, Blen, ctx->minfo);
        flint_free(A->exps);
        A->exps = newAexps;
        A->bits = newAbits;
    }

    _fmpz_vec_scalar_mul_fmpz(A->coeffs, B->coeffs, Blen, &Ccoeff0);
    _fmpz_mpoly_set_length(A, Blen, ctx);

cleanup_C:

    if (freeCcoeff0)
        fmpz_clear(&Ccoeff0);
}

