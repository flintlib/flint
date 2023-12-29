/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"
#include "fq_zech_mpoly.h"

void fq_zech_mpoly_ctx_init_deg(fq_zech_mpoly_ctx_t ctx, slong nvars,
                                  const ordering_t ord, mp_limb_t p, slong deg)
{
    fmpz_t P;
    mpoly_ctx_init(ctx->minfo, nvars, ord);
    fmpz_init_set_ui(P, p);
    fq_zech_ctx_init(ctx->fqctx, P, deg, "#");
    fmpz_clear(P);
}
