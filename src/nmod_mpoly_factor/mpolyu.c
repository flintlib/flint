/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

int nmod_mpolyu_is_canonical(
    const nmod_mpolyu_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->length; i++)
    {
        if (!nmod_mpoly_is_canonical(A->coeffs + i, ctx))
            return 0;
        if (nmod_mpoly_is_zero(A->coeffs + i, ctx))
            return 0;
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}
