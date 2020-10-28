/*
    Copyright (C) 2018-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_sub(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                              const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t t1, t2;
    slong easy_length = B->zpoly->length + C->zpoly->length;

    if (fmpq_mpoly_is_zero(B, ctx))
    {
        fmpq_mpoly_neg(A, C, ctx);
        return;
    }

    if (fmpq_mpoly_is_zero(C, ctx))
    {
        fmpq_mpoly_set(A, B, ctx);
        return;
    }

    fmpz_init(t1);
    fmpz_init(t2);

    fmpq_gcd_cofactors(A->content, t1, t2, B->content, C->content);

    fmpz_neg(t2, t2);
    fmpz_mpoly_scalar_fmma(A->zpoly, B->zpoly, t1, C->zpoly, t2, ctx->zctx);

    fmpz_clear(t1);
    fmpz_clear(t2);

    fmpq_mpoly_reduce_easy(A, easy_length, ctx);
}
