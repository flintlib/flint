/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_divrem_ideal(nmod_mpoly_struct ** Q, nmod_mpoly_t R,
               const nmod_mpoly_t A, nmod_mpoly_struct * const * B, slong len,
                                                    const nmod_mpoly_ctx_t ctx)
{
    /* TODO !!! */
    nmod_mpoly_divrem_ideal_monagan_pearce(Q, R, A, B, len, ctx);
}
