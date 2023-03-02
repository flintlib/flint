/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/* leave vec1 undefined but valid if division is not exact */
int _fmpz_vec_scalar_divides_fmpz(fmpz * vec1, const fmpz * vec2,
                                                    slong len2, const fmpz_t x)
{
    slong i;
    fmpz_t r;

    fmpz_init(r);

    for (i = 0; i < len2; i++)
    {
        fmpz_fdiv_qr(vec1 + i, r, vec2 + i, x);
        if (!fmpz_is_zero(r))
        {
            fmpz_clear(r);
            return 0;
        }
    }

    fmpz_clear(r);
    return 1;
}


int fmpz_mpoly_scalar_divides_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                    const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    int divides;
    slong N;

    if (A != B)
    {
        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        fmpz_mpoly_fit_length(A, B->length, ctx);
        fmpz_mpoly_fit_bits(A, B->bits, ctx);
        A->bits = B->bits;
        mpn_copyi(A->exps, B->exps, N*B->length);
    }
    divides = _fmpz_vec_scalar_divides_fmpz(A->coeffs, B->coeffs, B->length, c);
    _fmpz_mpoly_set_length(A, divides ? B->length : WORD(0), ctx);
    return divides;
}

int fmpz_mpoly_scalar_divides_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           ulong c, const fmpz_mpoly_ctx_t ctx)
{
    int divides;
    fmpz_t t;
    fmpz_init_set_ui(t, c);
    divides = fmpz_mpoly_scalar_divides_fmpz(A, B, t, ctx);
    fmpz_clear(t);
    return divides;
}

int fmpz_mpoly_scalar_divides_si(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           slong c, const fmpz_mpoly_ctx_t ctx)
{
    int divides;
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_si(t, c);
    divides = fmpz_mpoly_scalar_divides_fmpz(A, B, t, ctx);
    fmpz_clear(t);
    return divides;
}
