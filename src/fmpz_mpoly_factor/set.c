/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


void fmpz_mpoly_factor_set(
    fmpz_mpoly_factor_t A,
    const fmpz_mpoly_factor_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (A == B)
        return;

    fmpz_mpoly_factor_fit_length(A, B->num, ctx);
    fmpz_set(A->constant, B->constant);
    for (i = 0; i < B->num; i++)
    {
        fmpz_mpoly_set(A->poly + i, B->poly + i, ctx);
        fmpz_set(A->exp + i, B->exp + i);
    }
    A->num = B->num;
}
