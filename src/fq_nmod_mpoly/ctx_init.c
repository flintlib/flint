/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "mpoly.h"
#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_ctx_init_deg(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                                  const ordering_t ord, mp_limb_t p, slong deg)
{
    mpoly_ctx_init(ctx->minfo, nvars, ord);
    fq_nmod_ctx_init_ui(ctx->fqctx, p, deg, "#");
}

void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                               const ordering_t ord, const fq_nmod_ctx_t fqctx)
{
    mpoly_ctx_init(ctx->minfo, nvars, ord);
    fq_nmod_ctx_init_modulus(ctx->fqctx, fqctx->modulus, fqctx->var);
}
