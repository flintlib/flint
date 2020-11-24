/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_set_term_coeff_ui(nmod_mpoly_t A, slong i, ulong c,
                                                    const nmod_mpoly_ctx_t ctx)
{
    if (i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "nmod_mpoly_set_term_coeff_ui: index out of range.");
    }

    if (c >= ctx->mod.n)
    {
        NMOD_RED(c, c, ctx->mod);
    }

    A->coeffs[i] = c;
}
