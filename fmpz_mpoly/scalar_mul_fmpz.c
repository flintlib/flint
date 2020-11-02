/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                    const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if (fmpz_is_zero(c))
    {
        _fmpz_mpoly_set_length(A, 0, ctx);
        return;
    }

    if (A != B)
    {
        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        fmpz_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
        mpoly_copy_monomials(A->exps, B->exps, B->length, N);
    }

    if ((A != B) || !fmpz_is_one(c))
        _fmpz_vec_scalar_mul_fmpz(A->coeffs, B->coeffs, B->length, c);

    _fmpz_mpoly_set_length(A, B->length, ctx);
}

void fmpz_mpoly_scalar_mul_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           ulong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_ui(t, c);
    fmpz_mpoly_scalar_mul_fmpz(A, B, t, ctx);
    fmpz_clear(t);
}

void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           slong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_si(t, c);
    fmpz_mpoly_scalar_mul_fmpz(A, B, t, ctx);
    fmpz_clear(t);
}
