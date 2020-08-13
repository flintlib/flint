/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_scalar_mul_fmpz(fmpq_mpoly_t A,
              const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
{
    fmpq_mul_fmpz(A->content, B->content, c);
    if (fmpq_is_zero(A->content))
    {
        fmpz_mpoly_zero(A->zpoly, ctx->zctx);
    }
    else
    {
        fmpz_mpoly_set(A->zpoly, B->zpoly, ctx->zctx);
    }
}

void fmpq_mpoly_scalar_mul_ui(fmpq_mpoly_t A,
                     const fmpq_mpoly_t B, ulong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init_set_ui(C, c);
    fmpq_mpoly_scalar_mul_fmpz(A, B, C, ctx);
    fmpz_clear(C);
}

void fmpq_mpoly_scalar_mul_si(fmpq_mpoly_t A,
                     const fmpq_mpoly_t B, slong c, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t C;
    fmpz_init(C);
    fmpz_set_si(C, c);
    fmpq_mpoly_scalar_mul_fmpz(A, B, C, ctx);
    fmpz_clear(C);
}

