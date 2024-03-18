/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"
#include "mpoly.h"
#include "fq_zech_mpoly.h"

void fq_zech_mpoly_ctx_init_deg(fq_zech_mpoly_ctx_t ctx, slong nvars,
                                  const ordering_t ord, mp_limb_t p, slong deg)
{
    mpoly_ctx_init(ctx->minfo, nvars, ord);
    fq_zech_ctx_init_ui(ctx->fqctx, p, deg, "#");
}

void fq_zech_mpoly_ctx_clear(fq_zech_mpoly_ctx_t ctx)
{
    mpoly_ctx_clear(ctx->minfo);
    fq_zech_ctx_clear(ctx->fqctx);
}

void fq_zech_mpoly_ctx_change_modulus(fq_zech_mpoly_ctx_t ctx, slong deg)
{
    ulong p;
    p = fq_zech_ctx_mod(ctx->fqctx).n;
    fq_zech_ctx_clear(ctx->fqctx);
    fq_zech_ctx_init_ui(ctx->fqctx, p, deg, "#");
}
