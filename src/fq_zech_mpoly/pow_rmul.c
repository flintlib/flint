/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

void fq_zech_mpoly_pow_rmul(fq_zech_mpoly_t A, const fq_zech_mpoly_t B,
                                        ulong k, const fq_zech_mpoly_ctx_t ctx)
{
    fq_zech_mpoly_t T;
    fq_zech_mpoly_init(T, ctx);

    if (A == B)
    {
        fq_zech_mpoly_pow_rmul(T, A, k, ctx);
        fq_zech_mpoly_swap(T, A, ctx);
        goto cleanup;
    }

    fq_zech_mpoly_one(A, ctx);
    while (k > 0)
    { 
        fq_zech_mpoly_mul(T, A, B, ctx);
        fq_zech_mpoly_swap(A, T, ctx);
        k -= 1;
    }

cleanup:
    fq_zech_mpoly_clear(T, ctx);
}
