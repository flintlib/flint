/*
    Copyright (C) 2018-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_sub_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                    const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t t1, t2;
    slong easy_length = B->zpoly->length + 1;

    if (fmpq_is_zero(c))
    {
        fmpq_mpoly_set(A, B, ctx);
        return;
    }
    else if (fmpq_mpoly_is_zero(B, ctx))
    {
        fmpq_mpoly_set_fmpq(A, c, ctx);
        fmpq_mpoly_neg(A, A, ctx);
        return;
    }

    fmpz_init(t1);
    fmpz_init(t2);

    fmpq_gcd_cofactors(A->content, t1, t2, B->content, c);

    fmpz_mpoly_scalar_mul_fmpz(A->zpoly, B->zpoly, t1, ctx->zctx);
    fmpz_mpoly_sub_fmpz(A->zpoly, A->zpoly, t2, ctx->zctx);

    fmpz_clear(t1);
    fmpz_clear(t2);

    fmpq_mpoly_reduce_easy(A, easy_length, ctx);
}

void fmpq_mpoly_sub_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                    const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t t;
    *fmpq_numref(t) = *c;
    *fmpq_denref(t) = 1;
    fmpq_mpoly_sub_fmpq(A, B, t, ctx);
}

void fmpq_mpoly_sub_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                           ulong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t t;
    fmpz_init_set_ui(fmpq_numref(t), c);
    *fmpq_denref(t) = 1;
    fmpq_mpoly_sub_fmpq(A, B, t, ctx);
    fmpz_clear(fmpq_numref(t));
}

void fmpq_mpoly_sub_si(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                           slong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_t t;
    fmpz_init_set_si(fmpq_numref(t), c);
    *fmpq_denref(t) = 1;
    fmpq_mpoly_sub_fmpq(A, B, t, ctx);
    fmpz_clear(fmpq_numref(t));
}
