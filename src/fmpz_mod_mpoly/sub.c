/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

slong _fmpz_mod_mpoly_sub(
    fmpz * Acoeffs, ulong * Aexps,
    const fmpz * Bcoeffs, const ulong * Bexps, slong Blen,
    const fmpz * Ccoeffs, const ulong * Cexps, slong Clen,
    slong N,
    const ulong * cmpmask,
    const fmpz_mod_ctx_t fctx)
{
    slong i = 0, j = 0, k = 0;

    while (i < Blen && j < Clen)
    {
        int cmp = mpoly_monomial_cmp(Bexps + i*N, Cexps + j*N, N, cmpmask);

        if (cmp > 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            fmpz_set(Acoeffs + k, Bcoeffs + i);
            i++;
            k++;
        }
        else if (cmp == 0)
        {
            mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
            fmpz_mod_sub(Acoeffs + k, Bcoeffs + i, Ccoeffs + j, fctx);
            k += !fmpz_is_zero(Acoeffs + k);
            i++;
            j++;
        }
        else
        {
            mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
            fmpz_sub(Acoeffs + k, fmpz_mod_ctx_modulus(fctx), Ccoeffs + j);
            j++;
            k++;
        }
    }

    while (i < Blen)
    {
        mpoly_monomial_set(Aexps + k*N, Bexps + i*N, N);
        fmpz_set(Acoeffs + k, Bcoeffs + i);
        i++;
        k++;
    }

    while (j < Clen)
    {
        mpoly_monomial_set(Aexps + k*N, Cexps + j*N, N);
        fmpz_sub(Acoeffs + k, fmpz_mod_ctx_modulus(fctx), Ccoeffs + j);
        j++;
        k++;
    }

    return k;
}

void fmpz_mod_mpoly_sub(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_t C,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong Abits, N;
    ulong * Bexps = B->exps, * Cexps = C->exps;
    ulong * cmpmask;
    int freeBexps = 0, freeCexps = 0;
    TMP_INIT;

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        fmpz_mod_mpoly_neg(A, C, ctx);
        return;
    }
    else if (fmpz_mod_mpoly_is_zero(C, ctx))
    {
        fmpz_mod_mpoly_set(A, B, ctx);
        return;
    }

    TMP_START;
    Abits = FLINT_MAX(B->bits, C->bits);
    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    if (Abits != B->bits)
    {
        freeBexps = 1;
        Bexps = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexps, Abits, B->exps, B->bits,
                                                    B->length, ctx->minfo);
    }

    if (Abits != C->bits)
    {
        freeCexps = 1;
        Cexps = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexps, Abits, C->exps, C->bits,
                                                    C->length, ctx->minfo);
    }

    if (A == B || A == C)
    {
        fmpz_mod_mpoly_t T;
        fmpz_mod_mpoly_init3(T, B->length + C->length, Abits, ctx);
        T->length = _fmpz_mod_mpoly_sub(T->coeffs, T->exps, 
                                        B->coeffs, Bexps, B->length,
                                        C->coeffs, Cexps, C->length,
                                                      N, cmpmask, ctx->ffinfo);
        fmpz_mod_mpoly_swap(A, T, ctx);
        fmpz_mod_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mod_mpoly_fit_length_reset_bits(A, B->length + C->length, Abits, ctx);
        A->length = _fmpz_mod_mpoly_sub(A->coeffs, A->exps, 
                                        B->coeffs, Bexps, B->length,
                                        C->coeffs, Cexps, C->length,
                                                      N, cmpmask, ctx->ffinfo);
    }

    if (freeBexps)
        flint_free(Bexps);

    if (freeCexps)
        flint_free(Cexps);

    TMP_END;
}
