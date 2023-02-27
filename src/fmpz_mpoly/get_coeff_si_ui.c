/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

slong fmpz_mpoly_get_coeff_si_ui(const fmpz_mpoly_t A,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    index = mpoly_monomial_index_ui(A->exps, A->bits, A->length,
                                                              exp, ctx->minfo);
    if (index < 0)
    {
        return 0;
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        return fmpz_get_si(A->coeffs + index);
    }
}
