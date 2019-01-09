/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, slong nvars,
                                                        mp_limb_t p, slong deg)
{
    fmpz_t P;

    mpoly_ctx_init(ctx->minfo, nvars, ORD_LEX);

    fmpz_init_set_ui(P, p);
    fq_nmod_ctx_init(ctx->fqctx, P, deg, "#");
    fmpz_clear(P);
}
