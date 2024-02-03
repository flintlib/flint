/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_ctx_change_modulus(fq_nmod_mpoly_ctx_t ctx, slong deg)
{
    ulong p = ctx->fqctx->mod.n;
    fq_nmod_ctx_clear(ctx->fqctx);
    fq_nmod_ctx_init_ui(ctx->fqctx, p, deg, "#");
}
