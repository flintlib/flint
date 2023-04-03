/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/*
    c is assumed to be invertible and reduced mod n
*/
void nmod_mpoly_scalar_mul_nmod_invertible(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    mp_limb_t c,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(c != 0);
    FLINT_ASSERT(c < ctx->mod.n);
    FLINT_ASSERT(n_gcd(c, ctx->mod.n) == 1);

    if (A == B)
    {
        if (c == 1)
            return;
    }
    else
    {
        slong N;

        nmod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
        A->length = B->length;

        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        flint_mpn_copyi(A->exps, B->exps, N*B->length);
        if (c == UWORD(1))
        {
            flint_mpn_copyi(A->coeffs, B->coeffs, B->length);
            return;
        }
    }

    _nmod_vec_scalar_mul_nmod(A->coeffs, B->coeffs, B->length, c, ctx->mod);
}


/*
    c is assumed to be reduced mod n
*/
void nmod_mpoly_scalar_mul_nmod_general(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    mp_limb_t c,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, N;
    slong Alen, Blen;
    ulong * Aexp, * Bexp;
    mp_limb_t * Acoeff, * Bcoeff;

    FLINT_ASSERT(c < ctx->mod.n);

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    nmod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);

    Aexp = A->exps;
    Bexp = B->exps;
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;
    Blen = B->length;

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        mpoly_monomial_set(Aexp + N*Alen, Bexp + N*i, N);
        Acoeff[Alen] = nmod_mul(Bcoeff[i], c, ctx->mod);
        Alen += (Acoeff[Alen] != UWORD(0));
    }

    A->length = Alen;
}


void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                           ulong c, const nmod_mpoly_ctx_t ctx)
{
    if (c >= ctx->mod.n)
        NMOD_RED(c, c, ctx->mod);

    if (c == 0 || nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_zero(A, ctx);
        return;
    }

    if (n_gcd(c, ctx->mod.n) == 1)
        nmod_mpoly_scalar_mul_nmod_invertible(A, B, c, ctx);
    else
        nmod_mpoly_scalar_mul_nmod_general(A, B, c, ctx);

    return;
}
