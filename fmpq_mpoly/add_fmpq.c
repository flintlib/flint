/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_add_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                    const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
{

    if (fmpq_is_zero(c))
    {
        fmpq_mpoly_set(A, B, ctx);
        return;
    }
    else if (fmpq_mpoly_is_zero(B, ctx))
    {
        fmpq_mpoly_set_fmpq(A, c, ctx);
        return;
    }
    else
    {
        /*
            if B->content != 0 and c != 0,
                compute num/den = c / B->content
                A->content = B->content / den
                A->zpoly = den * B->zpoly + num
        */
        fmpq_t t;
        slong Blen = B->zpoly->length;

        fmpq_init(t);
        fmpq_div(t, c, B->content);
        if (!fmpz_is_one(fmpq_denref(t)))
        {
            fmpq_div_fmpz(A->content, B->content, fmpq_denref(t));
            fmpz_mpoly_scalar_mul_fmpz(A->zpoly, B->zpoly,
                                                    fmpq_denref(t), ctx->zctx);
            fmpz_mpoly_add_fmpz(A->zpoly, A->zpoly, fmpq_numref(t), ctx->zctx);
        }
        else
        {
            fmpq_set(A->content, B->content);
            fmpz_mpoly_add_fmpz(A->zpoly, B->zpoly, fmpq_numref(t), ctx->zctx);
        }
        /*
            optimization: since gcd(num, den) = 1 and B->zpoly is primitive,
                we do not need to reduce if the addition added a term
        */
        if (A->zpoly->length <= Blen)
        {
            fmpq_mpoly_reduce(A, ctx);
        }
        fmpq_clear(t);
        return;
    }
}

void fmpq_mpoly_add_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                    const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t t;
    fmpq_init(t);
    fmpz_set(fmpq_numref(t), c);
    fmpz_one(fmpq_denref(t));
    fmpq_mpoly_add_fmpq(A, B, t, ctx);
    fmpq_clear(t);
}

void fmpq_mpoly_add_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                           ulong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t t;
    fmpq_init(t);
    fmpz_set_ui(fmpq_numref(t), c);
    fmpz_one(fmpq_denref(t));
    fmpq_mpoly_add_fmpq(A, B, t, ctx);
    fmpq_clear(t);
}

void fmpq_mpoly_add_si(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                           slong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t t;
    fmpq_init(t);
    fmpz_set_si(fmpq_numref(t), c);
    fmpz_one(fmpq_denref(t));
    fmpq_mpoly_add_fmpq(A, B, t, ctx);
    fmpq_clear(t);
}
