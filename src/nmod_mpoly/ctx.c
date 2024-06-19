/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "mpoly.h"
#include "nmod_mpoly.h"

void nmod_mpoly_ctx_init(nmod_mpoly_ctx_t ctx, slong nvars,
                                       const ordering_t ord, ulong modulus)
{
    mpoly_ctx_init(ctx->minfo, nvars, ord);
    nmod_init(&ctx->mod, modulus);
}

void nmod_mpoly_ctx_init_rand(nmod_mpoly_ctx_t ctx, flint_rand_t state,
                                            slong max_nvars, ulong modulus)
{
    mpoly_ctx_init_rand(ctx->minfo, state, max_nvars);
    nmod_init(&ctx->mod, modulus);
}

void nmod_mpoly_ctx_set_modulus(nmod_mpoly_ctx_t ctx, ulong modulus)
{
    nmod_init(&ctx->mod, modulus);
}

void nmod_mpoly_ctx_clear(nmod_mpoly_ctx_t ctx)
{
    mpoly_ctx_clear(ctx->minfo);
}
