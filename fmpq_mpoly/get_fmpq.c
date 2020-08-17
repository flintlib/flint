/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_get_fmpq(fmpq_t c, const fmpq_mpoly_t A,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    /* the easiest way to throw is to call fmpz_mpoly_get_fmpz */
    fmpz_one(fmpq_denref(c));
    fmpz_mpoly_get_fmpz(fmpq_numref(c), A->zpoly, ctx->zctx); /* either 0 or 1 */
    fmpq_mul(c, c, A->content);
}
