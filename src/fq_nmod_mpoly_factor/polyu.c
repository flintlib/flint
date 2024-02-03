/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "n_poly.h"
#include "fq_nmod_mpoly_factor.h"

int n_polyu_fq_is_canonical(
    const n_polyu_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!n_fq_is_canonical(A->coeffs + d*i, ctx))
            return 0;
        if (_n_fq_is_zero(A->coeffs + d*i, d))
            return 0;
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}
