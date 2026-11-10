/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_mpoly.h"
#include "fmpq_mpoly.h"

int fmpq_mpoly_resultant(fmpq_mpoly_t R, const fmpq_mpoly_t A,
           const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)
{
    int success;
    fmpz_t degA, degB;

    /* the degrees are read before R, which may alias A or B, is written */
    fmpz_init(degA);
    fmpz_init(degB);
    fmpz_mpoly_degree_fmpz(degA, A->zpoly, var, ctx->zctx);
    fmpz_mpoly_degree_fmpz(degB, B->zpoly, var, ctx->zctx);

    /* res(cA A, cB B) = cA^deg(B) cB^deg(A) res(A, B) */
    success = fmpz_mpoly_resultant(R->zpoly, A->zpoly, B->zpoly, var, ctx->zctx);

    if (success && fmpz_sgn(degA) >= 0 && fmpz_sgn(degB) >= 0)
    {
        fmpq_t t;
        fmpq_init(t);

        success = fmpq_pow_fmpz(t, A->content, degB) &&
                  fmpq_pow_fmpz(R->content, B->content, degA);
        if (success)
            fmpq_mul(R->content, R->content, t);

        fmpq_clear(t);
    }
    else
    {
        fmpq_zero(R->content);
    }

    fmpq_mpoly_reduce(R, ctx);

    fmpz_clear(degA);
    fmpz_clear(degB);

    return success;
}
