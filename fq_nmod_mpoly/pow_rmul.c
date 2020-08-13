/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_pow_rmul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                        ulong k, const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpoly_t T;
    fq_nmod_mpoly_init(T, ctx);

    if (A == B)
    {
        fq_nmod_mpoly_pow_rmul(T, A, k, ctx);
        fq_nmod_mpoly_swap(T, A, ctx);
        goto cleanup;
    }

    fq_nmod_mpoly_set_ui(A, 1, ctx);
    while (k >= 1)
    { 
        fq_nmod_mpoly_mul_johnson(T, A, B, ctx);
        fq_nmod_mpoly_swap(A, T, ctx);
        k -= 1;
    }

cleanup:
    fq_nmod_mpoly_clear(T, ctx);
}
