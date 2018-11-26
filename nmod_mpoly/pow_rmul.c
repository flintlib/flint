/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_pow_rmul(nmod_mpoly_t A, const nmod_mpoly_t B,
                                         ulong k, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_t T;
    nmod_mpoly_init(T, ctx);

    if (A == B)
    {
        nmod_mpoly_pow_rmul(T, A, k, ctx);
        nmod_mpoly_swap(T, A, ctx);
        goto cleanup;
    }

    nmod_mpoly_set_ui(A, 1, ctx);
    while (k >= 1)
    { 
        nmod_mpoly_mul_johnson(T, A, B, ctx);
        nmod_mpoly_swap(A, T, ctx);
        k -= 1;
    }

cleanup:
    nmod_mpoly_clear(T, ctx);
}
