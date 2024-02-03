/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fq_nmod.h"

void fq_nmod_init(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_poly_init2_preinv(rop, ctx->mod.n, ctx->mod.ninv, fq_nmod_ctx_degree(ctx));
}

void fq_nmod_init2(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_poly_init2_preinv(rop, ctx->mod.n, ctx->mod.ninv, fq_nmod_ctx_degree(ctx));
}

void fq_nmod_clear(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
{
    nmod_poly_clear(rop);
}
