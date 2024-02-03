/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "mpoly.h"
#include "fmpz_mod_mpoly.h"

int fmpz_mod_mpoly_is_gen(const fmpz_mod_mpoly_t A,
                                     slong var, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (A->length != 1)
        return fmpz_is_one(fmpz_mod_ctx_modulus(ctx->ffinfo));

    if (!fmpz_is_one(A->coeffs + 0))
        return 0;

    return mpoly_is_gen(A->exps, var, A->bits, ctx->minfo);
}
