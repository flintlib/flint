/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_nmod_mpoly_factor.h"

int fq_nmod_mpoly_factor_cmp(
    const fq_nmod_mpoly_factor_t A,
    const fq_nmod_mpoly_factor_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    int cmp;
    slong i;

    cmp = fq_nmod_cmp(A->constant, B->constant, ctx->fqctx);
    if (cmp != 0)
        return cmp;

    if (A->num != B->num)
        return A->num > B->num ? 1 : -1;

    for (i = 0; i < A->num; i++)
    {
        cmp = fmpz_cmp(A->exp + i, B->exp + i);
        if (cmp != 0)
            return cmp;

        cmp = fq_nmod_mpoly_cmp(A->poly + i, B->poly + i, ctx);
        if (cmp != 0)
            return cmp;
    }

    return 0;
}
