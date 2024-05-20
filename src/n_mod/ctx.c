/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_mod.h"

void n_mod_ctx_init_rand(n_mod_ctx_ptr ctx, flint_rand_t state)
{
    ulong n = N_MOD_MIN + n_randint(state, N_MOD_MAX - N_MOD_MIN + 1);
    unsigned int norm = flint_clz(n);

    ctx->nu = n;
    ctx->nn = n << norm;
    ctx->ninv = n_preinvert_limb_prenorm(ctx->nn);
    ctx->norm = norm;
}
