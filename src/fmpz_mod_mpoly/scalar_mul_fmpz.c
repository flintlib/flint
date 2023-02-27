/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"
#include "fmpz_mod_mpoly.h"

/* c is assumed to be reduced and invertible mod n */
void fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N, Blen = B->length;

    FLINT_ASSERT(fmpz_mod_is_invertible(c, ctx->ffinfo));

    if (A != B)
    {
        fmpz_mod_mpoly_fit_length_reset_bits(A, Blen, B->bits, ctx);
        A->length = Blen;

        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        mpoly_copy_monomials(A->exps, B->exps, Blen, N);
    }

    _fmpz_mod_vec_scalar_mul_fmpz_mod(A->coeffs, B->coeffs, Blen, c, ctx->ffinfo);
}

/* c is assumed to be reduced mod n */
void fmpz_mod_mpoly_scalar_mul_fmpz_mod(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, N;
    slong Alen;
    slong Blen = B->length;
    ulong * Aexp;
    const ulong * Bexp;
    fmpz * Acoeff;
    const fmpz * Bcoeff;

    if (fmpz_is_zero(c) || Blen < 1)
    {
        fmpz_mod_mpoly_zero(A, ctx);
        return;
    }

    if (fmpz_is_one(c) || (Blen > 10 && fmpz_mod_is_invertible(c, ctx->ffinfo)))
    {
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(A, B, c, ctx);
        return;
    }

    N = mpoly_words_per_exp(B->bits, ctx->minfo);

    fmpz_mod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);

    Aexp = A->exps;
    Bexp = B->exps;
    Acoeff = A->coeffs;
    Bcoeff = B->coeffs;

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        mpoly_monomial_set(Aexp + N*Alen, Bexp + N*i, N);
        fmpz_mod_mul(Acoeff + Alen, Bcoeff + i, c, ctx->ffinfo);
        Alen += !fmpz_is_zero(Acoeff + Alen);
    }

    A->length = Alen;
}


void fmpz_mod_mpoly_scalar_mul_fmpz(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_is_canonical(c, ctx->ffinfo))
    {
        fmpz_mod_mpoly_scalar_mul_fmpz_mod(A, B, c, ctx);
    }
    else
    {
        fmpz_t cc;
        fmpz_init(cc);
        fmpz_mod_set_fmpz(cc, c, ctx->ffinfo);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod(A, B, cc, ctx);
        fmpz_clear(cc);
    }
}

void fmpz_mod_mpoly_scalar_mul_ui(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    ulong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t cc;
    fmpz_init(cc);
    fmpz_mod_set_ui(cc, c, ctx->ffinfo);
    fmpz_mod_mpoly_scalar_mul_fmpz_mod(A, B, cc, ctx);
    fmpz_clear(cc);
}

void fmpz_mod_mpoly_scalar_mul_si(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    slong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t cc;
    fmpz_init(cc);
    fmpz_mod_set_si(cc, c, ctx->ffinfo);
    fmpz_mod_mpoly_scalar_mul_fmpz_mod(A, B, cc, ctx);
    fmpz_clear(cc);
}

