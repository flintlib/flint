/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

int fmpz_mpoly_is_gen(const fmpz_mpoly_t A,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    if (A->length != WORD(1))
        return 0;

    if (!fmpz_is_one(A->coeffs + 0))
        return 0;

    return mpoly_is_gen(A->exps, var, A->bits, ctx->minfo);
}
