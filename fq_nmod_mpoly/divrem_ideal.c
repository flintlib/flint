/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_divrem_ideal(fq_nmod_mpoly_struct ** Q, fq_nmod_mpoly_t R,
         const fq_nmod_mpoly_t A, fq_nmod_mpoly_struct * const * B, slong len,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    /* nothing fancy */
    fq_nmod_mpoly_divrem_ideal_monagan_pearce(Q, R, A, B, len, ctx);
}
