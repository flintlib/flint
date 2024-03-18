/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "fq_nmod_mpoly_factor.h"

void fq_nmod_mpoly_factor_one(fq_nmod_mpoly_factor_t a, const fq_nmod_mpoly_ctx_t ctx)
{
	fq_nmod_one(a->constant, ctx->fqctx);
	a->num = 0;
}
