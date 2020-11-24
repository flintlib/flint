/*
    Copyright (C) 2017, 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_make_monic(nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    if (B->length == 0)
    {
        flint_throw(FLINT_ERROR, "nmod_mpoly_make_monic: polynomial is zero.");
    }

    nmod_mpoly_scalar_mul_nmod_invertible(A, B,
                                        nmod_inv(B->coeffs[0], ctx->mod), ctx);
}
