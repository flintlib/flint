/*
    Copyright (C) 2018, 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t A,
                                fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    index = mpoly_monomial_index_pfmpz(A->exps, A->bits, A->length,
                                                              exp, ctx->minfo);
    if (index < 0)
    {
        fmpz_zero(c);
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        fmpz_set(c, A->coeffs + index);
    }
}
