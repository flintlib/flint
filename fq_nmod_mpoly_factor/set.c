/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"

void fq_nmod_mpoly_factor_set(
    fq_nmod_mpoly_factor_t a,
    const fq_nmod_mpoly_factor_t b,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (a == b)
        return;

    fq_nmod_mpoly_factor_fit_length(a, b->num, ctx);
    fq_nmod_set(a->constant, b->constant, ctx->fqctx);
    for (i = 0; i < b->num; i++)
    {
        fq_nmod_mpoly_set(a->poly + i, b->poly + i, ctx);
        fmpz_set(a->exp + i, b->exp + i);
    }
    a->num = b->num;
}
