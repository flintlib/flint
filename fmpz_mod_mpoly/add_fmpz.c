/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_add_fmpz_mod(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong Blen = B->length;

    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx->ffinfo));

    if (fmpz_is_zero(c))
    {
        fmpz_mod_mpoly_set(A, B, ctx);
        return;
    }

    if (fmpz_mod_mpoly_is_zero(B, ctx))
    {
        fmpz_mod_mpoly_set_fmpz_mod(A, c, ctx);
        return;
    }

    if (mpoly_monomial_is_zero(B->exps + (Blen - 1)*N, N))
    {
        if (A != B)
        {
            fmpz_mod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
            _fmpz_vec_set(A->coeffs, B->coeffs, Blen - 1);
            mpoly_copy_monomials(A->exps, B->exps, Blen, N);
            _fmpz_mod_mpoly_set_length(A, B->length, ctx);
        }

        fmpz_mod_add(A->coeffs + Blen - 1, B->coeffs + Blen - 1, c, ctx->ffinfo);
        if (fmpz_is_zero(A->coeffs + Blen - 1))
            _fmpz_mod_mpoly_set_length(A, Blen - 1, ctx);
    }
    else
    {
        if (A != B)
        {
            fmpz_mod_mpoly_fit_length_reset_bits(A, Blen + 1, B->bits, ctx);
            _fmpz_vec_set(A->coeffs, B->coeffs, Blen);
            mpoly_copy_monomials(A->exps, B->exps, Blen, N);
        }
        else
        {
            fmpz_mod_mpoly_fit_length(A, Blen + 1, ctx);
        }

        mpoly_monomial_zero(A->exps + N*Blen, N);
        fmpz_set(A->coeffs + Blen, c);
        _fmpz_mod_mpoly_set_length(A, Blen + 1, ctx);
    }
}


void fmpz_mod_mpoly_add_fmpz(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t cc;

    if (fmpz_mod_is_canonical(c, ctx->ffinfo))
    {
        fmpz_mod_mpoly_add_fmpz_mod(A, B, c, ctx);
        return;
    }

    fmpz_init(cc);
    fmpz_mod_set_fmpz(cc, c, ctx->ffinfo);
    fmpz_mod_mpoly_add_fmpz_mod(A, B, cc, ctx);
    fmpz_clear(cc);
}

void fmpz_mod_mpoly_add_ui(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    ulong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t cc;
    fmpz_init(cc);
    fmpz_mod_set_ui(cc, c, ctx->ffinfo);
    fmpz_mod_mpoly_add_fmpz_mod(A, B, cc, ctx);
    fmpz_clear(cc);
}

void fmpz_mod_mpoly_add_si(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    slong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t cc;
    fmpz_init(cc);
    fmpz_mod_set_si(cc, c, ctx->ffinfo);
    fmpz_mod_mpoly_add_fmpz_mod(A, B, cc, ctx);
    fmpz_clear(cc);
}
