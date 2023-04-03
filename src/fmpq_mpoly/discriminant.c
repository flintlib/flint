/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_discriminant(fmpq_mpoly_t R, const fmpq_mpoly_t A,
                                     slong var, const fmpq_mpoly_ctx_t ctx)
{
    int success;
    fmpz_mpoly_univar_t Ax;

    fmpz_mpoly_univar_init(Ax, ctx->zctx);

    fmpz_mpoly_to_univar(Ax, A->zpoly, var, ctx->zctx);

    success = fmpz_mpoly_univar_discriminant(R->zpoly, Ax, ctx->zctx);

    if (success && Ax->length > 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mul_ui(t, Ax->exps + 0, 2);
        fmpz_sub_ui(t, t, 2);
        success = fmpq_pow_fmpz(R->content, A->content, t);
        fmpz_clear(t);
    }
    else
    {
        fmpq_zero(R->content);
    }

    fmpq_mpoly_reduce(R, ctx);

    fmpz_mpoly_univar_clear(Ax, ctx->zctx);

    return success;
}

