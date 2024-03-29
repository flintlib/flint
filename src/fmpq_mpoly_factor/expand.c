/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly_factor.h"


int fmpq_mpoly_factor_expand(
    fmpq_mpoly_t A,
    const fmpq_mpoly_factor_t f,
    const fmpq_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i;
    fmpq_mpoly_t t1, t2;

    fmpq_mpoly_init(t1, ctx);
    fmpq_mpoly_init(t2, ctx);

    fmpq_mpoly_set_fmpq(A, f->constant, ctx);

    for (i = 0; i < f->num; i++)
    {
        if (fmpz_sgn(f->exp + i) < 0 ||
            !fmpq_mpoly_pow_fmpz(t1, f->poly + i, f->exp + i, ctx))
        {
            success = 0;
            goto cleanup;
        }
        fmpq_mpoly_mul(t2, A, t1, ctx);
        fmpq_mpoly_swap(A, t2, ctx);
    }

cleanup:

    fmpq_mpoly_clear(t1, ctx);
    fmpq_mpoly_clear(t2, ctx);

    return success;
}
