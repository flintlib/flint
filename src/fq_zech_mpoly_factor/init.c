/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly_factor.h"

void fq_zech_mpoly_factor_init(
    fq_zech_mpoly_factor_t f,
    const fq_zech_mpoly_ctx_t ctx)
{
    fq_zech_init(f->constant, ctx->fqctx);
	fq_zech_one(f->constant, ctx->fqctx);
    f->poly  = NULL;
    f->exp   = NULL;
    f->num   = 0;
    f->alloc = 0;
}
