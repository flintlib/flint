/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_get_fmpq(fmpq_t x, const fmpq_mpoly_t poly,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    fmpz_one(fmpq_denref(x));
    fmpz_mpoly_get_fmpz(fmpq_numref(x), poly->zpoly, ctx->zctx); /* either 0 or 1 */
    fmpq_mul(x, x, poly->content);
}
