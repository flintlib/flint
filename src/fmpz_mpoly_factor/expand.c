/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mpoly_factor.h"

int fmpz_mpoly_factor_expand(
    fmpz_mpoly_t A,
    const fmpz_mpoly_factor_t f,
    const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i;
    fmpz_mpoly_t t1, t2;

    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);

    fmpz_mpoly_set_fmpz(A, f->constant, ctx);

    for (i = 0; i < f->num; i++)
    {
        if (fmpz_sgn(f->exp + i) < 0 ||
            !fmpz_mpoly_pow_fmpz(t1, f->poly + i, f->exp + i, ctx))
        {
            success = 0;
            goto cleanup;
        }
        fmpz_mpoly_mul(t2, A, t1, ctx);
        fmpz_mpoly_swap(A, t2, ctx);
    }

cleanup:

    fmpz_mpoly_clear(t1, ctx);
    fmpz_mpoly_clear(t2, ctx);

    return success;
}
