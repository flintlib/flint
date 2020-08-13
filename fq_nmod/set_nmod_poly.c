/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_set_nmod_poly(fq_nmod_t a, const nmod_poly_t b,
                                                       const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(a->mod.n == b->mod.n);
    FLINT_ASSERT(a->mod.n == ctx->modulus->mod.n);

    if (b->length <= 2*(ctx->modulus->length - 1))
    {
        nmod_poly_set(a, b);
        fq_nmod_reduce(a, ctx);
    }
    else
    {
        nmod_poly_rem(a, b, ctx->modulus);
    }
}
