/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


int fmpz_mpoly_factor_add(
    fmpz_mpoly_factor_t A,
    const fmpz_mpoly_factor_t B,
    const fmpz_mpoly_factor_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    int success = 0;
    fmpz_mpoly_t b, c;

    fmpz_mpoly_init(b, ctx);
    fmpz_mpoly_init(c, ctx);

    if (!fmpz_mpoly_factor_expand(b, B, ctx))
        goto cleanup;
    if (!fmpz_mpoly_factor_expand(c, C, ctx))
        goto cleanup;

    fmpz_mpoly_factor_fit_length(A, 1, ctx);
    fmpz_one(A->constant);
    fmpz_mpoly_add(A->poly + 0, b, c, ctx);
    fmpz_one(A->exp + 0);
    A->num = 1;
    success = 1;

cleanup:

    fmpz_mpoly_clear(c, ctx);
    fmpz_mpoly_clear(b, ctx);
    return success;
}
