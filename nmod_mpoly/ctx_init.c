/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_ctx_init(nmod_mpoly_ctx_t ctx, slong nvars,
                                           const ordering_t ord, ulong modulus)
{
    ctx->n = (ord == ORD_DEGLEX || ord == ORD_DEGREVLEX) ? nvars + 1 : nvars;
    ctx->ord = ord;
    nmodf_ctx_init(ctx->ffinfo, modulus);
}
