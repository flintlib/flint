/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_pow_rmul(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B,
                                        ulong k, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_t T;
    fmpz_mod_mpoly_init(T, ctx);

    if (A == B)
    {
        fmpz_mod_mpoly_pow_rmul(T, A, k, ctx);
        fmpz_mod_mpoly_swap(T, A, ctx);
        goto cleanup;
    }

    fmpz_mod_mpoly_one(A, ctx);
    while (k >= 1)
    { 
        fmpz_mod_mpoly_mul(T, A, B, ctx);
        fmpz_mod_mpoly_swap(A, T, ctx);
        k -= 1;
    }

cleanup:
    fmpz_mod_mpoly_clear(T, ctx);
}
