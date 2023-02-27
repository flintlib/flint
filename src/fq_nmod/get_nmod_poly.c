/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_get_nmod_poly(nmod_poly_t a, const fq_nmod_t b,
                                                       const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(b->mod.n == ctx->modulus->mod.n);
    a->mod = ctx->modulus->mod;
    nmod_poly_set(a, b);
}
