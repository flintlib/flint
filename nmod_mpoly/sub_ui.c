/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_sub_ui(nmod_mpoly_t A, const nmod_mpoly_t B,
                                           ulong c, const nmod_mpoly_ctx_t ctx)
{
    if (c >= ctx->mod.n)
    {
        NMOD_RED(c, c, ctx->mod);
    }
    nmod_mpoly_add_ui(A, B, nmod_neg(c, ctx->mod), ctx);
}
