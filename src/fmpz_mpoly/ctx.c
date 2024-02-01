/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, slong nvars, const ordering_t ord)
{
    mpoly_ctx_init(ctx->minfo, nvars, ord);
}

void fmpz_mpoly_ctx_init_rand(fmpz_mpoly_ctx_t ctx, flint_rand_t state, slong max_nvars)
{
    mpoly_ctx_init_rand(ctx->minfo, state, max_nvars);
}

void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx)
{
    mpoly_ctx_clear(ctx->minfo);
}
