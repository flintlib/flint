/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_is_fmpq(fmpq_t x, const fmpq_mpoly_t poly,
                                                   const fmpq_mpoly_ctx_t ctx)
{
    if (fmpq_mpoly_is_zero(poly, ctx)) {
        fmpq_zero(x);
        return 1;
    } else if (fmpz_mpoly_is_one(poly->zpoly, ctx->zctx))
    {
        fmpq_set(x, poly->content);
        return 1;
    } else
    {
        return 0;
    }
}
