/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_resultant(fmpq_mpoly_t R, const fmpq_mpoly_t A,
           const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpoly_univar_t Ax, Bx;

    fmpz_mpoly_univar_init(Ax, ctx->zctx);
    fmpz_mpoly_univar_init(Bx, ctx->zctx);

    fmpz_mpoly_to_univar(Ax, A->zpoly, var, ctx->zctx);
    fmpz_mpoly_to_univar(Bx, B->zpoly, var, ctx->zctx);

    success = fmpz_mpoly_univar_resultant(R->zpoly, Ax, Bx, ctx->zctx);

    if (success && Ax->length > 0 && Bx->length > 0)
    {
        fmpq_t t;
        fmpq_init(t);

        success = fmpq_pow_fmpz(t, A->content, Bx->exps + 0) &&
                  fmpq_pow_fmpz(R->content, B->content, Ax->exps + 0);
        if (success)
            fmpq_mul(R->content, R->content, t);

        fmpq_clear(t);
    }
    else
    {
        fmpq_zero(R->content);
    }

    fmpq_mpoly_reduce(R, ctx);

    fmpz_mpoly_univar_clear(Ax, ctx->zctx);
    fmpz_mpoly_univar_clear(Bx, ctx->zctx);

    return success;
}

