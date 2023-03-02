/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

void nmod_mpoly_factor_set(nmod_mpoly_factor_t res,
                     const nmod_mpoly_factor_t fac, const nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (res == fac)
        return;

    nmod_mpoly_factor_fit_length(res, fac->num, ctx);
    res->constant = fac->constant;
    for (i = 0; i < fac->num; i++)
    {
        nmod_mpoly_set(res->poly + i, fac->poly + i, ctx);
        fmpz_set(res->exp + i, fac->exp + i);
    }
    res->num = fac->num;
}
