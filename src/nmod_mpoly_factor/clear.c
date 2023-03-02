/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"

void nmod_mpoly_factor_clear(nmod_mpoly_factor_t fac, const nmod_mpoly_ctx_t ctx)
{
    if (fac->alloc > 0)
    {
        slong i;

        for (i = 0; i < fac->alloc; i++)
        {
            nmod_mpoly_clear(fac->poly + i, ctx);
			fmpz_clear(fac->exp + i);
        }

        flint_free(fac->poly);
        flint_free(fac->exp);
    }
}
