/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017, 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_add_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                    const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
{
    slong i, N;
    slong Blen = B->length;

    if (Blen == 0)
    {
        fmpz_mpoly_set_fmpz(A, c, ctx);
        return;
    }

    if (!fmpz_is_zero(c))
    {
        N = mpoly_words_per_exp(B->bits, ctx->minfo);

        if (mpoly_monomial_is_zero(B->exps + (Blen - 1)*N, N))
        {
            if (A != B)
            {
                fmpz_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);

                for (i = 0; i < Blen - 1; i++)
                    fmpz_set(A->coeffs + i, B->coeffs + i);

                for (i = 0; i < Blen*N; i++)
                    A->exps[i] = B->exps[i];

                _fmpz_mpoly_set_length(A, B->length, ctx);
            }

            fmpz_add(A->coeffs + Blen - 1, B->coeffs + Blen - 1, c);

            if (fmpz_is_zero(A->coeffs + Blen - 1))
                _fmpz_mpoly_set_length(A, Blen - 1, ctx);
        }
        else
        {
            fmpz_mpoly_fit_length(A, Blen + 1, ctx);

            if (A != B)
            {
                fmpz_mpoly_fit_bits(A, B->bits, ctx);
                A->bits = B->bits;

                for (i = 0; i < Blen; i++)
                    fmpz_set(A->coeffs + i, B->coeffs + i);

                for (i = 0; i < Blen*N; i++)
                    A->exps[i] = B->exps[i];
            } 

            for (i = 0; i < N; i++)
                A->exps[Blen*N + i] = 0;

            fmpz_set(A->coeffs + Blen, c);

            _fmpz_mpoly_set_length(A, Blen + 1, ctx);
        }
    }
    else if (A != B)
    {
        fmpz_mpoly_set(A, B, ctx);
    }
}

void fmpz_mpoly_add_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           ulong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_ui(t, c);
    fmpz_mpoly_add_fmpz(A, B, t, ctx);
    fmpz_clear(t);
}

void fmpz_mpoly_add_si(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                           slong c, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_si(t, c);
    fmpz_mpoly_add_fmpz(A, B, t, ctx);
    fmpz_clear(t);
}
