/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_mul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                        const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)
{
    /* nothing fancy for now */
    fq_nmod_mpoly_mul_johnson(A, B, C, ctx);
}
