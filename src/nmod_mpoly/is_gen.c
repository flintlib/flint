/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int nmod_mpoly_is_gen(const nmod_mpoly_t A,
                                         slong var, const nmod_mpoly_ctx_t ctx)
{
    if (A->length != WORD(1))
        return 0;

    if (A->coeffs[0] != UWORD(1))
        return 0;

    return mpoly_is_gen(A->exps, var, A->bits, ctx->minfo);
}
