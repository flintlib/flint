/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fq.h"

void _fq_sparse_reduce(fmpz *R, slong lenR, const fq_ctx_t ctx)
{
    const slong d = ctx->j[ctx->len - 1];

    FMPZ_VEC_NORM(R, lenR);

    if (lenR > d)
    {
        slong i, k;

        for (i = lenR - 1; i >= d; i--)
        {
            for (k = ctx->len - 2; k >= 0; k--)
            {
                fmpz_submul(R + ctx->j[k] + i - d, R + i, ctx->a + k);
            }
            fmpz_zero(R + i);
        }
    }

    _fmpz_vec_scalar_mod_fmpz(R, R, FLINT_MIN(d, lenR), fq_ctx_prime(ctx));
}
