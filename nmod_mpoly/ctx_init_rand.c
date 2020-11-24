/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_ctx_init_rand(nmod_mpoly_ctx_t ctx, flint_rand_t state,
                                            slong max_nvars, mp_limb_t modulus)
{
    mpoly_ctx_init_rand(ctx->minfo, state, max_nvars);
    nmod_init(&ctx->mod, modulus);
}
